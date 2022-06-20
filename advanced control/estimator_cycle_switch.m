function x_estim = estimator_cycle_switch(process_time,cycle_time,...
               stations_working,u,measurements,operating_vars,x_estim,...
               n_cycle,control_mode, sampling_interval, control_interval)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function estimating parameters and states at every carousel rotation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs
    %
    % process_time      =   timer started at process onset (s)
    % cycle_time        =   timer re-started at every carousel rotation (s)
    % stations_working     =   vector [1x4] - for i=1:4:
    %                       - stations_working(i)=1 if station i was processing material during cycle that just finished ;
    %                       - stations_working(i)=0 if station i was empty during cycle that just finished; 
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
    %                       Fields follow the structure defined in run_carousel.m
    % n_cycle           =   cycle counter - number of cycle that has just finished
    % control_mode      =   scalar defined in run_carousel.m
    % sampling_interval     =   sampling interval (s)
    % control_interval  =   control_interval (s)
    %
    % Outputs           
    % x_estim           =   Updated after call to this function 
    %                       Object containing states and parameters
    %                       estimated by estimator_online.m and estimator_cycle_switch
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    if control_mode>=3.2
        if n_cycle > 0
            if stations_working(1)==1
                %% Parameter estimation filtration section 
                % estimated cake length, slurry volume and slurry concentration
                % imposed equal to respective measurements
                x_estim.cake_counter=x_estim.cake_counter+1;
                t_exp=measurements.t_meas(end-round(u.t_cycle/sampling_interval):end);
                t_exp=t_exp-t_exp(1);
                x_estim.(['batch_' num2str(n_cycle)]).cake_number=x_estim.cake_counter;
                x_estim.(['batch_' num2str(n_cycle)]).L_cake=...
                    measurements.L_cake_LI101(end);
                x_estim.(['batch_' num2str(n_cycle)]).V_slurry=...
                    measurements.V_slurry_LI101(end);
                x_estim.(['batch_' num2str(n_cycle)]).c_slurry=...
                    measurements.c_slurry_AI101(end);

                % filtrate mass profile retrieval
                m_exp=measurements.m_filt_WI101(end-round(u.t_cycle/sampling_interval):end);
                V_exp=(m_exp-m_exp(1))/842;

                % estimation of porosity
                A=181.4584e-6;
                x_estim.(['batch_' num2str(n_cycle)]).E=1-x_estim.(['batch_' num2str(n_cycle)]).V_slurry*...
                    x_estim.(['batch_' num2str(n_cycle)]).c_slurry/1293/...
                    x_estim.(['batch_' num2str(n_cycle)]).L_cake/A;

                % calculate approximate filtration duration
                m_solid_initial=x_estim.(['batch_' num2str(n_cycle)]).V_slurry*...
                    x_estim.(['batch_' num2str(n_cycle)]).c_slurry;
                V_filt_final=x_estim.(['batch_' num2str(n_cycle)]).V_slurry*(1-...
                    x_estim.(['batch_' num2str(n_cycle)]).c_slurry/1293)-m_solid_initial/...
                    1293*x_estim.(['batch_' num2str(n_cycle)]).E/(1-x_estim.(['batch_' num2str(n_cycle)]).E);
                c=(m_solid_initial)/V_filt_final; % Mass of dry cake deposited per unit volume of filtrate [kg sol/ m3 filtrate liquid]more compact way, but same results as Brigi's
                filt_duration=(1.1339e-003*2.7e9*c*V_filt_final^2/(2*A^2*u.P_compr)+...
                    1.1339e-003*3e9*V_filt_final/(A*u.P_compr));     

                % estimation of specific cake resistance and filter mesh resistance
                options=optimoptions('fmincon','Algorithm', 'sqp', 'StepTolerance', 1e-9,...
                    'OptimalityTolerance',1e-9,'Display','off');
                x0=[2.7e9 5e9]/1e9;
                [opt,~,~,~,~,~,hessian]  = fmincon(@obj, x0,...
                    [],[],[],[],[0.1 0.1],[20 20],[],options,x_estim, u, t_exp, V_exp, n_cycle, filt_duration);            
                % estimated alpha
                x_estim.(['batch_' num2str(n_cycle)]).alpha(1)=opt(1)*1e9;
                % estimated filter mesh resistance that ['batch_' num2str(n_cycle)]
                % encounters in Stations 1-4 (current resistance in
                % Stations 2-4 will be different, but when ['batch_' num2str(n_cycle)]
                % arrives there, this estimation should be alright 
                x_estim.(['batch_' num2str(n_cycle)]).Rm=opt(2)*1e9;
                % store estimation uncertainty
                err = sqrt(diag(inv(hessian)))*1e9;
                x_estim.(['batch_' num2str(n_cycle)]).err_alpha=err(1);
                x_estim.(['batch_' num2str(n_cycle)]).err_Rm=err(2);
                % updated filtration duration based on estimated specific cake
                % resistance and filter mesh resistance
                x_estim.(['batch_' num2str(n_cycle)]).filt_duration=(1.1339e-003*...
                    x_estim.(['batch_' num2str(n_cycle)]).alpha*...
                    c*V_filt_final^2/(2*A^2*u.P_compr)+1.1339e-003*...
                    x_estim.(['batch_' num2str(n_cycle)]).Rm(1)*V_filt_final/(A*u.P_compr)); 
            else
                x_estim.(['batch_' num2str(n_cycle)])='no loaded material';
            end
            %% Parameter estimation drying section 
            % Estimation of initial conditions before of drying
            if stations_working(3)==1                               
                duration_deliq=sum(operating_vars.t_cycle_vector(end-1:end))+u.t_cycle-...
                    x_estim.(['batch_' num2str(n_cycle-2)]).filt_duration;             
                deliq_output = layer2_model_deliquoring(x_estim.(['batch_' num2str(n_cycle-2)]),duration_deliq,u);
                x_estim.(['batch_' num2str(n_cycle-2)]).S_initial_drying=deliq_output.S;
                x_estim.(['batch_' num2str(n_cycle-2)]).t_estim=0;     
                
            %% Initialization of EKF
                if control_mode>=3.2
                    p_ekf.number_nodes=10;  
                    p_ekf.step_grid_drying=x_estim.(['batch_' num2str(n_cycle-2)]).L_cake/p_ekf.number_nodes;
                    p_ekf.nodes_drying=linspace(p_ekf.step_grid_drying/2, ...
                        x_estim.(['batch_' num2str(n_cycle-2)]).L_cake-p_ekf.step_grid_drying/2,...
                        p_ekf.number_nodes);                
                    p_ekf.station_diameter = 0.0152;
                    p_ekf.A = p_ekf.station_diameter.^2.*pi./4;     
                    p_ekf.MW_components=46.07*1e-3;
                    p_ekf.MW_air = 28.971e-3;    
                    cp_air = [29e-3, 0.2199e-5, 0.5723e-8, -2.871e-12]/p_ekf.MW_air*1e3; % nitrogen specific heat [J/(kg K)] 
                    cp_EtOH= [61.34e-3, 15.72e-5, 8.749e-8, 19.83e-12]/p_ekf.MW_components*1e3;
                    p_ekf.cp_gas_components=[cp_air; cp_EtOH];
                    p_ekf.visc_gas_phase = 1.8e-5;          % air viscosity at room temperature and 1 atm [Pa s] (temperature effect considered in drying model)
                    p_ekf.cp=@(T) cp_air(1)+cp_air(2)*(T+273.15)+cp_air(3)*(T+273.15)^2+cp_air(4)*(T+273.15)^3;
                    p_ekf.a_V = 126000;  
                    p_ekf.zeta=199.78;
                    p_ekf.rho_sol = 1293;                % PCM density [kg/m^3]    
                    p_ekf.cp_s = 2067; 
                    p_ekf.vl_crit = 0.05;                    % Drying - critical impurity content [m3_i/m3]
                    p_ekf.vl_eq = 0.0005;                       % Drying - equilibrium impurity content [m3_l/m3] 
                    p_ekf.control_interval=control_interval;
                    p_ekf.drying_sampling_interval=sampling_interval;
                    p_ekf.coeff_antoine= [74.475, -7164.3,-7.327, 3.1340e-6, 2];
                    p_ekf.rho_liq_components=842;
                    p_ekf.latent_heat = 846*1e3;  
                    p_ekf.cp_liq_components = 2570;                       
                    p_ekf.L_cake=x_estim.(['batch_' num2str(n_cycle-2)]).L_cake;
                    p_ekf.E=x_estim.(['batch_' num2str(n_cycle-2)]).E;
                    p_ekf.rho_liq_components=842;
                    p_ekf.rho_sol=1293;
                    p_ekf.T_room=295.25;

                    Tg_0=ones(1,p_ekf.number_nodes)*295.25;
                    Ts_0=ones(1,p_ekf.number_nodes)*295.25;
                    mass_frG_0=zeros(1,p_ekf.number_nodes);                               

                    S0=x_estim.(['batch_' num2str(n_cycle-2)]).S_initial_drying;
                    epsL_0=S0*x_estim.(['batch_' num2str(n_cycle-2)]).E;
                    vl_0=epsL_0;  

                    % Filter construction - scaled ekf
                    x0=[Tg_0'/300; Ts_0'/300; mass_frG_0'/0.01; reshape(vl_0',...
                        size(vl_0,1)*size(vl_0,2),1)/0.01]';                
                    S1=triu(tril(ones(p_ekf.number_nodes)),-1); % sparsity matrix - consistent with the FVM upwind differencing scheme lower diagonal        
                    S=repmat(S1,[4,4]);
                    p_ekf.options = odeset('JPattern',S);                  
                    filter = extendedKalmanFilter(@ekf_ModelDrying_scaled,...
                        @ekf_MeasurementFcn,x0); 
                    filter.StateCovariance = diag([zeros(1,3*p_ekf.number_nodes) ones(1,p_ekf.number_nodes)*5e-3]);
                    filter.MeasurementNoise = eye(1)*3e-6; 
                    filter.ProcessNoise = diag([ones(1,p_ekf.number_nodes) ones(1,p_ekf.number_nodes) ones(1,p_ekf.number_nodes) ones(1,p_ekf.number_nodes)])*50e-6;

                    % estimated states vector
                    initial_epsl=x_estim.(['batch_' num2str(n_cycle-2)]).E*...
                        x_estim.(['batch_' num2str(n_cycle-2)]).S_initial_drying;
                    rho_cake=(initial_epsl*p_ekf.rho_liq_components+(1-...
                        x_estim.(['batch_' num2str(n_cycle-2)]).E)*p_ekf.rho_sol);
                    x_estim.(['batch_' num2str(n_cycle-2)]).mass_avg=mean(initial_epsl./...
                        rho_cake.*p_ekf.rho_liq_components);                
                    x_estim.(['batch_' num2str(n_cycle-2)]).mass_avg_std=sqrt(filter.StateCovariance(end))*0.01;          
                    x_estim.(['batch_' num2str(n_cycle-2)]).Tg_bot=Tg_0;
                    x_estim.(['batch_' num2str(n_cycle-2)]).Ts_bot=Ts_0;
                    x_estim.(['batch_' num2str(n_cycle-2)]).wg=mass_frG_0;
                    x_estim.(['batch_' num2str(n_cycle-2)]).Tg_top=p_ekf.T_room;
                    
                    % store EKF object
                    x_estim.(['batch_' num2str(n_cycle-2)]).EKF=filter;
                    x_estim.(['batch_' num2str(n_cycle-2)]).p=p_ekf;                    
                end
            end
        else
            x_estim.cake_counter=0;
        end
    end
end

function L=obj(x_opt, x_estim, u, t_exp, V_exp, n_rotation, filt_duration)
       
    [V_filt, V_exp]=filt_model(x_opt, x_estim, u, t_exp, V_exp, n_rotation,filt_duration);
    SSE=sum((V_filt-V_exp).^2);
    Nexp=length(V_filt);
    L=Nexp/2*log(SSE);
    
end

function [V_filt,V_exp]=filt_model(x_opt,x,u,t_exp,V_exp,n_cycle,filt_duration) 
    % variables
    alpha=x_opt(1)*1e9;
    Rm=x_opt(2)*1e9;

    % inputs and parameters
    A=181.4584e-6; % filter section (m2)
    visc_liq=1.1339e-003; % liquid viscosity (Pa s)
    dP=u.P_compr; %
    m_solid_initial=x.(['batch_' num2str(n_cycle)]).V_slurry*...
        x.(['batch_' num2str(n_cycle)]).c_slurry;
    V_filt_final=x.(['batch_' num2str(n_cycle)]).V_slurry*(1-...
        x.(['batch_' num2str(n_cycle)]).c_slurry/1293)-m_solid_initial/...
        1293*x.(['batch_' num2str(n_cycle)]).E/(1-x.(['batch_' num2str(n_cycle)]).E);
       
    % calculations
    c=(m_solid_initial)/V_filt_final; % Mass of dry cake deposited per unit volume of filtrate [kg sol/ m3 filtrate liquid]more compact way, but same results as Brigi's
    t_filt=t_exp;
    a_filt= visc_liq.*alpha.*c./(2.*A.^2.*dP);  % darcy coeff 1
    b_filt= visc_liq.*Rm./(A.*dP);      % darcy coeff 2
    V_filt= (-b_filt+sqrt(b_filt.^2+4.*a_filt.*t_filt))./(2.*a_filt);     % V(t) - filtrate volume [m^3]
    
    % approximate filtration duration not depending on x_opt, so that
    % number of points for SSE calculation does not change
    t_filt=t_filt(t_filt<=filt_duration);
    V_filt=V_filt(t_filt<=filt_duration);
    V_exp=V_exp(t_filt<=filt_duration);   
        
end

