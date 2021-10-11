function [x,y,measurements]=switch_cycle(~,cryst_output,p,d,u,x,y,measurements,n_cycle)
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

            % feed from position 1, cake properties, filtration and
            % deliquoring variables
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
            x.pos1.x_perc=cryst_output.x_perc; 
            x.pos1.CSD_perc=cryst_output.CSD_perc;
            
            % initial filtration states
            x.pos1.filtration_time=0;
            x.pos1.filtration_duration=1;
            x.pos1.V_filt=0;
            x.pos1.filtration_finished=0;

            % cake and liquid properties from slurry fed at current process time
            x.pos1.E=0.35.*d.E_dist;    % porosity_function_shape(p.station_diameter,x.pos1.CSD,x.pos1.x,1);           
%             x.pos1.m0=trapz(x.pos1.x,x.pos1.CSD); 
%             x.pos1.m1=trapz(x.pos1.x,x.pos1.CSD.*x.pos1.x); 
%             x.pos1.m2=trapz(x.pos1.x,x.pos1.CSD.*(x.pos1.x.^2)); 
%             x.pos1.m3=trapz(x.pos1.x,x.pos1.CSD.*(x.pos1.x.^3)); 
            x.pos1.alpha=2.7e9.*d.alpha_dist;   % sum(p.alpha_CSD(x.pos1.x_perc,x.pos1.E).*x.pos1.CSD_perc);
            x.pos1.k=1/(x.pos1.alpha*p.rho_sol*(1-x.pos1.E));    
            x.pos1.a_V=126000;  % 6*x.pos1.m2/x.pos1.m3;

            x.pos1.visc_liq=p.visc_liq_components(x.pos1.T);                        
            x.pos1.V_slurry_initial=u.V_slurry*d.V_slurry_dist*p.ports_working(1);        % Slurry volume fed at the beginning of the filtration batch [m^3]
            x.pos1.c_slurry=cryst_output.conc_MSMPR.*d.c_slurry_dist*p.ports_working(1);
            x.pos1.m_solid_initial=x.pos1.V_slurry_initial.*cryst_output.conc_MSMPR.*d.c_slurry_dist;    % Total solid mass filled into the filter [kg]
            x.pos1.V_solid_initial=x.pos1.m_solid_initial./p.rho_sol;        % Total solid volume filled into the filter [m^3]
            x.pos1.V_liq_initial=x.pos1.V_slurry_initial-x.pos1.V_solid_initial;  % Total liquid volume filled into the filter [m^3]
            x.pos1.V_liquid_pores_end_of_filtration=x.pos1.V_solid_initial*x.pos1.E/(1-x.pos1.E); % Volume of liquid in the pores at the end of filtration [m3]
            x.pos1.V_filt_final=x.pos1.V_liq_initial-x.pos1.V_liquid_pores_end_of_filtration;   % volume of filtrate at the end of filtration [m3]
            x.pos1.c=(x.pos1.m_solid_initial)/x.pos1.V_filt_final; % Mass of dry cake deposited per unit volume of filtrate [kg sol/ m3 filtrate liquid]more compact way, but same results as Brigi's
            x.pos1.L_cake=(x.pos1.m_solid_initial/p.rho_sol/(1-x.pos1.E))/p.A; % height of the deposited cake [m]    
            x.pos1.Pb=sum(p.pb_CSD(x.pos1.x_perc,x.pos1.E).*x.pos1.CSD_perc);
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
            y.pos1.(['batch_' num2str(n_cycle)]).P=1e5;
        end

        %% Sensors vectors creation
        if n_cycle == 1
            measurements.t=0;
            measurements.m_filt_WI101=0;
            measurements.P_PI102=1e5;
            measurements.c_slurry_AI101=250;
            measurements.L_cake_LI101=0;
            measurements.V_slurry_LI101=0;            
            measurements.Tg_top_TI101=u.Tinlet_drying;
            measurements.Tg_bot_TI102=295.2;
            measurements.Vdryer_FI101=0;
            
        end                
        x.m_filt_bias=measurements.m_filt_WI101(end);
        x.c_slurry_vector(n_cycle)=x.pos1.c_slurry;
end
