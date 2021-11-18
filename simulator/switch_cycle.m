function [x,y,measurements,measurements_nf,operating_vars]=switch_cycle(t,cryst_output,p,d,u,...
                x,y,measurements,measurements_nf,operating_vars,n_cycle)
        % This function handles the routines to be called at the end of each cycle:
        % -     the first time that some material enters a processing position (1-4), the
        %       measurement vectors for that position is created
        % -     the states in each position are transferred to the
        %       following one, e.g.: x.pos3 = x.pos2
        % -     the properties of the cake and of the liquid phase obtained 
        %       with the slurry entering the carousel the carousel at the
        %       current time step are calculated and become part of the
        %       state vector. Proceeding this way, disturbances in the feed
        %       occurring in the future will not affect batches that
        %       entered the carousel previously
        
        %% Add sensors readings and operating variable profiles during carousel rotation and mesh cleaning
        if t>0  
        % measurements
        sampling_times=measurements.t_meas(end):p.filtration_sampling_interval:t;
        measurements.t_meas=[measurements.t_meas sampling_times];
        measurements.m_filt_WI101=[measurements.m_filt_WI101 measurements.m_filt_WI101(end)*ones(1,length(sampling_times))+randn(1,length(sampling_times))*0.00005];
        measurements.P_PI101=[measurements.P_PI101 u.P_compr(end)*ones(1,length(sampling_times))];
        measurements.P_PI102=[measurements.P_PI102 0*ones(1,length(sampling_times))];
        measurements.c_slurry_AI101=[measurements.c_slurry_AI101 0*ones(1,length(sampling_times))];
        measurements.L_cake_LI101=[measurements.L_cake_LI101 0*ones(1,length(sampling_times))];
        measurements.V_slurry_LI101=[measurements.V_slurry_LI101 0*ones(1,length(sampling_times))];        
        measurements.Tg_in_TI101=[measurements.Tg_in_TI101 round(u.Tinlet_drying,1)*ones(1,length(sampling_times))];
        measurements.Tg_out_TI102=[measurements.Tg_out_TI102 round(p.T_room,1)*ones(1,length(sampling_times))];
        measurements.Vdryer_FI101=[measurements.Vdryer_FI101 0*ones(1,length(sampling_times))];
% 
        % noise-free measurements
        measurements_nf.t_meas=[measurements_nf.t_meas sampling_times];
        measurements_nf.m_filt_WI101=[measurements_nf.m_filt_WI101 measurements_nf.m_filt_WI101(end)*ones(1,length(sampling_times))];
        measurements_nf.P_PI101=[measurements_nf.P_PI101 u.P_compr(end)*ones(1,length(sampling_times))];
        measurements_nf.P_PI102=[measurements_nf.P_PI102 0*ones(1,length(sampling_times))];
        measurements_nf.c_slurry_AI101=[measurements_nf.c_slurry_AI101 0*ones(1,length(sampling_times))];
        measurements_nf.L_cake_LI101=[measurements_nf.L_cake_LI101 0*ones(1,length(sampling_times))];
        measurements_nf.V_slurry_LI101=[measurements_nf.V_slurry_LI101 0*ones(1,length(sampling_times))];        
        measurements_nf.Tg_in_TI101=[measurements_nf.Tg_in_TI101 round(u.Tinlet_drying,1)];
        measurements_nf.Tg_out_TI102=[measurements_nf.Tg_out_TI102 round(p.T_room,1)];
        measurements_nf.Vdryer_FI101=[measurements_nf.Vdryer_FI101 0];
        
        % operating variables
            sampling_times_control=operating_vars.t_vector(end):p.control_interval:t;
            operating_vars.t_vector=[operating_vars.t_vector sampling_times_control];     
            operating_vars.P_compr_vector=[operating_vars.P_compr_vector u.P_compr*ones(1,length(sampling_times_control))];
            operating_vars.Tin_drying_vector=[operating_vars.Tin_drying_vector u.Tinlet_drying*ones(1,length(sampling_times_control))];

        end
        %% Station 4
        if n_cycle > 3

            % measurements
            y.pos4.(['batch_' num2str(n_cycle-3)]).t=0;
            % if some deliquoring occurs
            y.pos4.(['batch_' num2str(n_cycle-3)]).m_filt=0;

            % drying
            y.pos4.(['batch_' num2str(n_cycle-3)]).Vdryer=0;            
            y.pos4.(['batch_' num2str(n_cycle-3)]).w_EtOH_gas=zeros(p.number_volatile_components,1);            
            initial_epsl=x.pos3.E*x.pos3.S;
            y.pos4.(['batch_' num2str(n_cycle-3)]).S=mean(x.pos3.S);
            rho_cake=(initial_epsl*p.rho_liq_components+(1-x.pos3.E)*p.rho_sol);
            y.pos4.(['batch_' num2str(n_cycle-3)]).w_EtOH_cake=mean(initial_epsl./rho_cake.*p.rho_liq_components);
            y.pos4.(['batch_' num2str(n_cycle-3)]).Tg_bot=x.pos3.T;
            y.pos4.(['batch_' num2str(n_cycle-3)]).Ts_bot=x.pos3.T;
            y.pos4.(['batch_' num2str(n_cycle-3)]).Tg_top=x.pos3.T;

            
            % feed from position 1, cake properties, filtration and
            % deliquoring variables
            x.pos4=x.pos3;
            
            % create new states for drying
            x.pos4.Tg=ones(1,x.pos4.number_nodes_drying)*x.pos4.T;
            x.pos4.Ts=x.pos4.Tg; 
            x.pos4.gas_mass_fr_vect=ones(1,p.number_volatile_components*x.pos4.number_nodes_drying)*0;
            x.pos4.Tin=x.pos4.T;                         
        end
        %% Station 3
        if n_cycle > 2

            % measurements
            y.pos3.(['batch_' num2str(n_cycle-2)]).t=0;
            y.pos3.(['batch_' num2str(n_cycle-2)]).m_filt=0;
            y.pos3.(['batch_' num2str(n_cycle-2)]).S=mean(x.pos2.S);
            initial_epsl=x.pos2.E*x.pos2.S;
            rho_cake=(initial_epsl*p.rho_liq_components+(1-x.pos2.E)*p.rho_sol);
            y.pos3.(['batch_' num2str(n_cycle-2)]).w_EtOH_cake=mean(initial_epsl./rho_cake.*p.rho_liq_components);

            x.pos3=x.pos2;
        end

        %% Station 2
        if n_cycle > 1

            % measurements
            y.pos2.(['batch_' num2str(n_cycle-1)]).t=0;
            y.pos2.(['batch_' num2str(n_cycle-1)]).m_filt=0;
            y.pos2.(['batch_' num2str(n_cycle-1)]).S=mean(x.pos1.S);
            initial_epsl=x.pos1.E*x.pos1.S;
            rho_cake=(initial_epsl*p.rho_liq_components+(1-x.pos1.E)*p.rho_sol);
            y.pos2.(['batch_' num2str(n_cycle-1)]).w_EtOH_cake=mean(initial_epsl./rho_cake.*p.rho_liq_components);   
  
            % feed from position 1, cake properties, filtration and
            % deliquoring variables
            x.pos2=x.pos1;           
        end
        
        %% Station 1
        if n_cycle > 0            
            % feed     
            x.pos1.T=cryst_output.T;
            x.pos1.x=cryst_output.x;%*u.CSD_dist;
            x.pos1.CSD=cryst_output.CSD;
            x.pos1.CSD_perc=cryst_output.CSD_perc;
            
            % initial filtration states
            x.pos1.filtration_time=0;
            x.pos1.filtration_duration=1;
            x.pos1.V_filt=0;
            x.pos1.filtration_finished=0;
            
            % cake and liquid properties from slurry fed at current process time
            x.pos1.E=0.35.*d.E_dist;    %            
            x.pos1.alpha=2.7e9.*d.alpha_dist;   % 
            x.pos1.k=1/(x.pos1.alpha*p.rho_sol*(1-x.pos1.E));    
            x.pos1.a_V=126000;  % 1/m

            x.pos1.visc_liq=p.visc_liq_components(x.pos1.T);                        
            x.pos1.V_slurry_initial=u.V_slurry*d.V_slurry_dist*p.stations_working(1);        % Slurry volume fed at the beginning of the filtration batch [m^3]
            x.pos1.c_slurry=cryst_output.conc_slurry.*d.c_slurry_dist*p.stations_working(1);
            x.pos1.m_solid_initial=x.pos1.V_slurry_initial.*cryst_output.conc_slurry.*d.c_slurry_dist;    % Total solid mass filled into the filter [kg]
            x.pos1.V_solid_initial=x.pos1.m_solid_initial./p.rho_sol;        % Total solid volume filled into the filter [m^3]
            x.pos1.V_liq_initial=x.pos1.V_slurry_initial-x.pos1.V_solid_initial;  % Total liquid volume filled into the filter [m^3]
            x.pos1.V_liquid_pores_end_of_filtration=x.pos1.V_solid_initial*x.pos1.E/(1-x.pos1.E); % Volume of liquid in the pores at the end of filtration [m3]
            x.pos1.V_filt_final=x.pos1.V_liq_initial-x.pos1.V_liquid_pores_end_of_filtration;   % volume of filtrate at the end of filtration [m3]
            x.pos1.c=(x.pos1.m_solid_initial)/x.pos1.V_filt_final; % Mass of dry cake deposited per unit volume of filtrate [kg sol/ m3 filtrate liquid]more compact way, but same results as Brigi's
            x.pos1.L_cake=(x.pos1.m_solid_initial/p.rho_sol/(1-x.pos1.E))/p.A; % height of the deposited cake [m]    
            x.pos1.Pb=sum(p.pb_CSD(x.pos1.x,x.pos1.E).*x.pos1.CSD_perc);
            x.pos1.V_liquid_pores=x.pos1.V_liquid_pores_end_of_filtration;
            x.pos1.number_nodes_deliq=max(round(x.pos1.L_cake/p.min_length_discr)+1,2);
            x.pos1.number_nodes_drying=x.pos1.number_nodes_deliq;
            x.pos1.step_grid_deliq=x.pos1.L_cake/x.pos1.number_nodes_deliq;
            x.pos1.nodes_deliq=linspace(x.pos1.step_grid_deliq/2,x.pos1.L_cake-x.pos1.step_grid_deliq/2,x.pos1.number_nodes_deliq);
            x.pos1.step_grid_drying=x.pos1.L_cake/x.pos1.number_nodes_drying;
            x.pos1.nodes_drying=linspace(x.pos1.step_grid_drying/2,...
                x.pos1.L_cake-x.pos1.step_grid_drying/2, x.pos1.number_nodes_drying);
            
            % outputs Station 1 measurements initialization
            y.pos1.(['batch_' num2str(n_cycle)]).t=0;
            y.pos1.(['batch_' num2str(n_cycle)]).m_filt=0;
            y.pos1.(['batch_' num2str(n_cycle)]).S=1;
            initial_epsl=x.pos1.E;
            rho_cake=(initial_epsl*p.rho_liq_components+(1-x.pos1.E)*p.rho_sol);
            y.pos1.(['batch_' num2str(n_cycle)]).w_EtOH_cake=mean(initial_epsl./rho_cake.*p.rho_liq_components);
            y.pos1.(['batch_' num2str(n_cycle)]).L_cake=x.pos1.L_cake;
            y.pos1.(['batch_' num2str(n_cycle)]).c_slurry=x.pos1.c_slurry;
            y.pos1.(['batch_' num2str(n_cycle)]).V_slurry=x.pos1.V_slurry_initial;
            % Processing operation sequence indicators
            y.sequence.(['batch_' num2str(n_cycle)]).filtration_duration=0;
            y.sequence.(['batch_' num2str(n_cycle)]).deliquoring_duration=0;
            y.sequence.(['batch_' num2str(n_cycle)]).drying_duration=0;
            
        end

          
        
        %% Initializations
        x.m_filt_bias=measurements.m_filt_WI101(end);
        x.c_slurry_vector(n_cycle)=x.pos1.c_slurry;
end
