function [x,y]=switch_phase(t,cryst_output,p,u,x,y,n_rotation)
        
        %% Position 5
        if n_rotation > 3
            % measurements
            % if some deliquoring occurs
            y.pos4.(['cycle_' num2str(n_rotation-3)]).t_filt=0;
            y.pos4.(['cycle_' num2str(n_rotation-3)]).V_filt=0;
            y.pos4.(['cycle_' num2str(n_rotation-3)]).w_filt=...
                y.pos3.(['cycle_' num2str(n_rotation-3)]).w_filt(:,end);
            
            % drying
            y.pos4.(['cycle_' num2str(n_rotation-3)]).t_drying=0;
            y.pos4.(['cycle_' num2str(n_rotation-3)]).wg=zeros(p.number_volatile_components,1);
            y.pos4.(['cycle_' num2str(n_rotation-3)]).Tg=298;
                
            % feed from position 1, cake properties, filtration and
            % deliquoring variables
            x.pos4=x.pos3;
            x.pos4.biasV =0; % correction for washing filtrate volume measurement vector

            
            % create new states for drying
            x.pos4.Tg=ones(1,x.pos4.number_nodes_drying)*x.pos4.T;
            x.pos4.Ts=x.pos4.Tg; 
            x.pos4.gas_mass_fr_vect=ones(1,p.number_volatile_components*x.pos4.number_nodes_drying)*0;
%             MW_gas=p.MW_N2;
%             for i = 1:p.number_volatile_components
%                 coeff=p.coeff_antoine(i,:);
%                 Psat_volatiles = (10.^(coeff(:,1)'-coeff(:,2)'./(coeff(:,3)'+(343-273.15)))*133.322); % Antoine Equation - saturation pressure [Pa]   
%                 mf(i)=Psat_volatiles*p.MW_components(i)/MW_gas./101325;
%                 
%             end            
%             MW_gas=mf*p.MW_components(1:2)+(1-sum(mf))*p.MW_N2;
%             mean(x.pos4.
        end


        %% Position 3
        if n_rotation > 2
            % measurements
            y.pos3.(['cycle_' num2str(n_rotation-2)]).t_filt=0;
            y.pos3.(['cycle_' num2str(n_rotation-2)]).V_filt=0;
            y.pos3.(['cycle_' num2str(n_rotation-2)]).w_filt=...
                y.pos2.(['cycle_' num2str(n_rotation-2)]).w_filt(:,end);

            % feed from position 1, cake properties, filtration and
            % deliquoring variables
            x.pos3=x.pos2;
            x.pos3.biasV =0; % correction for washing filtrate colume measurement vector

        end

        %% Position 2: washing
        if n_rotation > 1
            % measurements
            y.pos2.(['cycle_' num2str(n_rotation-1)]).t_filt=0;
            y.pos2.(['cycle_' num2str(n_rotation-1)]).V_filt=0;
            y.pos2.(['cycle_' num2str(n_rotation-1)]).w_filt=...
                y.pos1.(['cycle_' num2str(n_rotation-1)]).w_filt(:,end);

            % feed from position 1, cake properties, filtration and
            % deliquoring variables
            x.pos2=x.pos1;
            
            % washing states
            x.pos2.washing_time=0;
            x.pos2.washing_duration=1;
            x.pos2.V_wash=0;
                        
            x.pos2.V_wash_final= x.pos2.E*p.A* x.pos2.L_cake*u.W;
%             tw_sim=x.pos2.V_wash_final/Q_wash;
            x.pos2.biasV =0;
            x.pos2.washing_finished=0;

        end
        
        %% Position 1: filtration 
        if n_rotation > 0
            % measurements
            y.pos1.(['cycle_' num2str(n_rotation)]).t_filt=0;
            y.pos1.(['cycle_' num2str(n_rotation)]).V_filt=0;
            y.pos1.(['cycle_' num2str(n_rotation)]).w_filt=cryst_output.liq_mass_fr_vect;

            % feed     
            x.pos1.T=cryst_output.T;
            x.pos1.liq_mass_fr_vect=cryst_output.liq_mass_fr_vect;
            x.pos1.x=cryst_output.x;
            x.pos1.CSD=cryst_output.CSD;

            % cake properties
            x.pos1.E=porosity_function_shape(p.station_diameter,x.pos1.CSD,x.pos1.x,1);           % Ouchiyama model
            x.pos1.m0=trapz(x.pos1.x,x.pos1.CSD); 
            x.pos1.m1=trapz(x.pos1.x,x.pos1.CSD.*x.pos1.x); 
            x.pos1.m2=trapz(x.pos1.x,x.pos1.CSD.*(x.pos1.x.^2)); 
            x.pos1.m3=trapz(x.pos1.x,x.pos1.CSD.*(x.pos1.x.^3)); 
            x.pos1.alpha=trapz(x.pos1.x,p.alpha_CSD(x.pos1.x,x.pos1.E).*x.pos1.CSD/x.pos1.m0);
            x.pos1.k=1/(x.pos1.alpha*p.rho_sol*(1-x.pos1.E));      % cake permeability [?]
            x.pos1.a_V=6*x.pos1.m2/x.pos1.m3;

            % filtration
            x.pos1.visc_liq=p.visc_liquid_phase_from_mass_fr(x.pos1.T,x.pos1.liq_mass_fr_vect);
            x.pos1.V_slurry_initial=x.pos0.charge_cell_volume;        % Slurry volume fed at the beginning of the filtration batch [m^3]
            x.pos1.m_solid_initial=x.pos1.V_slurry_initial.*cryst_output.conc_MSMPR;    % Total solid mass filled into the filter [kg]
            x.pos1.V_solid_initial=x.pos1.m_solid_initial./p.rho_sol;        % Total solid volume filled into the filter [m^3]
            x.pos1.V_liq_initial=x.pos1.V_slurry_initial-x.pos1.V_solid_initial;  % Total liquid volume filled into the filter [m^3]
            x.pos1.V_liquid_pores_end_of_filtration=x.pos1.V_solid_initial*x.pos1.E/(1-x.pos1.E); % Volume of liquid in the pores at the end of filtration [m3]
            x.pos1.V_filt_final=x.pos1.V_liq_initial-x.pos1.V_liquid_pores_end_of_filtration;   % volume of filtrate at the end of filtration [m3]
            x.pos1.c=(x.pos1.m_solid_initial)/x.pos1.V_filt_final; % Mass of dry cake deposited per unit volume of filtrate [kg sol/ m3 filtrate liquid]more compact way, but same results as Brigi's
            x.pos1.L_cake=(x.pos1.c*x.pos1.V_filt_final/p.rho_sol/(1-x.pos1.E))/p.A; % height of the deposited cake [m]    
            x.pos1.Pb=trapz(x.pos1.x,p.pb_CSD(x.pos1.x,x.pos1.E).*x.pos1.CSD/x.pos1.m0);
            x.pos1.V_liquid_pores=x.pos1.V_liquid_pores_end_of_filtration;
            x.pos1.number_nodes_deliq=max(round(x.pos1.L_cake/p.min_length_discr)+1,2);
%             x.pos1.number_nodes_filt=1;
            x.pos1.number_nodes_drying=x.pos1.number_nodes_deliq;
            x.pos1.number_nodes_washing=max(x.pos1.number_nodes_deliq,p.number_nodes_washing);
            % deliquoring grid
            x.pos1.step_grid_deliq=x.pos1.L_cake/x.pos1.number_nodes_deliq;
            x.pos1.nodes_deliq=linspace(x.pos1.step_grid_deliq/2,x.pos1.L_cake-x.pos1.step_grid_deliq/2,x.pos1.number_nodes_deliq);
            % washing grid
            x.pos1.step_grid_washing=x.pos1.L_cake/x.pos1.number_nodes_washing;
            x.pos1.nodes_washing=linspace(x.pos1.step_grid_washing/2,x.pos1.L_cake-x.pos1.step_grid_washing/2,x.pos1.number_nodes_washing);  
            % drying grid
            x.pos1.step_grid_drying=x.pos1.L_cake/x.pos1.number_nodes_drying;
            x.pos1.nodes_drying=linspace(x.pos1.step_grid_drying/2,...
                x.pos1.L_cake-x.pos1.step_grid_drying/2, x.pos1.number_nodes_drying);
            
            % initial filtration states
            x.pos1.filtration_time=0;
            x.pos1.filtration_duration=1;
            x.pos1.V_filt=0;
            x.pos1.filtration_finished=0;
            
                        
            % filtration biases needed for building measurement vectors
            x.pos1.biasV_filt=0;
            x.pos1.biast_filt=0;

        end
        
        
        %% Position 0 - charge cell
        x.pos0.charge_cell_volume = 0;
        y.pos0.(['cycle_' num2str(n_rotation)]).V_filt = 0;
        
        if n_rotation ==1
            y.cont_sign.pos1_4.t=[0 t];
            y.cont_sign.pos1_4.V=[0 0];
        end


end
