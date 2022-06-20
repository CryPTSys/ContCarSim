function x_estim = estimator_online(process_time,cycle_time,...
               stations_working,u,measurements,operating_vars,x_estim,...
               n_cycle,control_mode, sampling_interval, control_interval,...
               simulation_step)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs
    %
    % process_time      =   timer started at process onset (s)
    % cycle_time        =   timer re-started at every carousel rotation (s)
    % stations_working     =   vector [1x4] - for i=1:4:
    %                       - stations_working(i)=1 if during current cycle station i is processing material;
    %                       - stations_working(i)=0 if during current cycle station i is empty; 
    % u                 =   vector of manipulated variables during previous control interval
    %                       Fields of u:
    %                       - u.t_cycle=cycle duration (s)
    %                       - u.V_slurry=fed slurry volume (m3)
    %                       - u.P_compr= gauge pressure provided by compressor P101 (Pa)
    %                       - u.Tinlet_drying=drying gas temperature Station 4 (K)                      
    % u_nominal         =   nominal value of manipulated variables, as set in run_carousel.m
    %                       Same fields of u
    % measurements      =   object of process measurements since process onset 
    %                       with the sampling interval that has been set in run_carousel.m
    %                       Fields of measurements:
    %                       - measurements.t_meas = vector of sampling intervals (s)
    %                       - measurements.m_filt_WI101 = vector of filtrate mass measured by WI101 (kg)
    %                       - measurements.P_PI101 = vector of pressure measured by PI101 - gauge (Pa)
    %                       - measurements.P_PI102 = vector of pressure measured by PI102 - gauge (Pa)
    %                       - measurements.c_slurry_AI101 = vector of slurry concentration measured by AI101 (kg/m3)
    %                       - measurements.L_cake_LI101 = vector of height of cakes in Station 1 measured by LI101 (m)
    %                       - measurements.V_slurry_LI101 = vector of slurry volume in Station 1 measured by LI101 (m)
    %                       - measurements.Tg_in_TI101 = vector of temperatures of drying gas measured by TI101 (K) - inlet
    %                       - measurements.Tg_out_TI102 = vector of temperatures of drying gas measured by TI102 (K) - outlet                      
    %                       - measurements.Vdryer_FI101 = vector of drying gas flowrate measured by FI101 (m3/s)
    % operating_vars  =   object storing the profiles of the manipulated variables (automatically updated)
    %                       Fields of operating_vars:
    %                       - operating_vars.t_vector = control times vector
    %                       - operating_vars.P_compr_vector = u.P_compr time profile [1 x length(operating_vars.t_vector)]
    %                       - operating_vars.Tin_drying_vector = u.Tinlet_drying time profile [1 x length(operating_vars.t_vector)]
    %                       - operating_vars.n_cycle_vector = list of number of completed carousel cycles
    %                       - operating_vars.t_cycle_vector = u.t_cycle time profile [1 x length(operating_vars.n_cycle_vector)]
    %                       - operating_vars.V_slurry_vector = u.V_slurry time profile [1 x length(operating_vars.n_cycle_vector)]
    % x_estim           =   object containing states and parameters estimated by estimator_online.m and estimator_cycle_switch
    % n_cycle           =   cycle counter - number of current cycle 
    % control_mode      =   scalar defined in run_carousel.m
    % sampling_interval     =   sampling interval (s)
    % control_interval  =   control interval (s)
    % simulation_step   =   duration of current simulation step (s)
    %
    % Outputs           
    % x_estim           =   Updated after call to this function 
    %                       Object containing states and parameters
    %                       estimated by estimator_online.m and estimator_cycle_switch
    %                       Fields follow the structure defined in run_carousel.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    % Run EKF
    if stations_working(4)==1 && control_mode>=3.2
        % retrieve filter object, input and output data
        filter=x_estim.(['batch_' num2str(n_cycle-3)]).EKF;
%         t_final=x_estim.(['batch_' num2str(n_cycle-3)]).t_estim(end)+...
%             simulation_step;
%         Tin=measurements.Tg_top_TI101(end-simulation_step/sampling_interval:end);
        Tout=measurements.Tg_out_TI102(end-simulation_step/sampling_interval:10:end);
        Vdryer=measurements.Vdryer_FI101(end-simulation_step/sampling_interval:10:end);
        time_steps_simulation=linspace(0, simulation_step, length(Tout));
        number_nodes=10;
        
        % update inlet temperature after HL
        heat_loss_pars=[0.6998,68.88485571,-0.7178];   
        a_hl=heat_loss_pars(1);
        tau=heat_loss_pars(2);
        ug_exp=heat_loss_pars(3);
        a=.155./(tau+.155);               
        time_steps_HL=linspace(0,simulation_step,max(ceil(simulation_step/.155),2));
        Tin(1)=x_estim.(['batch_' num2str(n_cycle-3)]).Tg_top(end);
        Tin_C=Tin(1)-273.15;
        ug0=Vdryer(2)*1e-3/60/x_estim.(['batch_' num2str(n_cycle-3)]).p.A;
        for i = 2:length(time_steps_HL)
            Tin_C=(1-a)*Tin_C+a*(u.Tinlet_drying-273.15-(a_hl+ug0*...
              ug_exp)*(u.Tinlet_drying-x_estim.(['batch_' num2str(n_cycle-3)]).p.T_room));              
            Tin(i)=Tin_C+273.15;
        end        
        Tin=interp1(time_steps_HL,Tin,time_steps_simulation);  
        
        % run EKF
        for i = 1:1:length(time_steps_simulation)-1   
            data.u=u;   
            data.Vdryer=Vdryer(i);
            data.Tin=Tin(i);
            data.p=x_estim.(['batch_' num2str(n_cycle-3)]).p;
            step=time_steps_simulation(i+1)-time_steps_simulation(i);

            % update process noise variance covariance matrix
            % (only if residual ethanol at cake top > equilibrium)
            if filter.State(31)*0.01>5e-3 
                filter.ProcessNoise = Qcalc(filter.State, data); 
            end
            
            % predict and update
            [PredStates,StCov] = predict(filter,data,step); % prediction step
            [CorrStates,StCov] = correct(filter, Tout(i+1)/300, 10);     
            vol_cont_impurities=CorrStates(31:40)';
            % correct numerical errors
            for k = 2:9
                if vol_cont_impurities(k)<vol_cont_impurities(k-1) || vol_cont_impurities(k)>vol_cont_impurities(k+1)
                    CorrStates(3*number_nodes+k)=PredStates(3*number_nodes+k);
                end
            end
            filter.State(filter.State(3*number_nodes+1:end)<5e-4)=5e-4;
            filter.State=CorrStates;

            % calculate and store
            x_estim.(['batch_' num2str(n_cycle-3)]).Tg_bot=[x_estim.(['batch_' num2str(n_cycle-3)]).Tg_bot; ...
                CorrStates(1:number_nodes)*300];
            x_estim.(['batch_' num2str(n_cycle-3)]).Ts_bot=[x_estim.(['batch_' num2str(n_cycle-3)]).Ts_bot; ...
                CorrStates(1+number_nodes:2*number_nodes)*300]; 
            x_estim.(['batch_' num2str(n_cycle-3)]).wg=[x_estim.(['batch_' num2str(n_cycle-3)]).wg; ...
                CorrStates(1+2*number_nodes:3*number_nodes)*0.01]; 
            vol_cont_impurities=CorrStates(1+3*number_nodes:4*number_nodes)'*0.01;  
            std_vol_cont_impurities=sqrt(diag(StCov(1+3*number_nodes:4*number_nodes,1+3*number_nodes:4*number_nodes)))*0.01;%(1+3*number_nodes:4*number_nodes))'*0.01;
            rho_cake=(vol_cont_impurities*842+(1-x_estim.(['batch_' num2str(n_cycle-3)]).E)*1293);
            x_estim.(['batch_' num2str(n_cycle-3)]).mass_avg=[x_estim.(['batch_' num2str(n_cycle-3)]).mass_avg,...
                mean(vol_cont_impurities./rho_cake.*842)];  
            x_estim.(['batch_' num2str(n_cycle-3)]).mass_avg_std=[x_estim.(['batch_' num2str(n_cycle-3)]).mass_avg_std,...
                min(mean(std_vol_cont_impurities./rho_cake*842),0.001)];
        end
        x_estim.(['batch_' num2str(n_cycle-3)]).t_estim=[x_estim.(['batch_' num2str(n_cycle-3)]).t_estim, ...
            x_estim.(['batch_' num2str(n_cycle-3)]).t_estim(end)+time_steps_simulation(2:end)];
        x_estim.(['batch_' num2str(n_cycle-3)]).Tg_top=[x_estim.(['batch_' num2str(n_cycle-3)]).Tg_top; ...
            Tin(2:end)'];
    end

end

