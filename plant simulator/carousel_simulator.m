function [x,y]=carousel_simulator(t,Dt,cryst_output,p,u,x,y,n_rotation)
        
        %% Slurry tank
%         x.slurry_tank_volume = x.slurry_tank_volume + ...
%             cryst_output.flowrate_MSMPR*Dt - u.flowrate_slurry*Dt;
    
        %% Position 0 - charge cell
        x.pos0.charge_cell_volume = u.V_slurry;%;x.pos0.charge_cell_volume + u.V_slurry;
        
        %% Position 1: filtration (+ pre-deliquoring_grad)
        if n_rotation > 0
            % update pressure drop through
            x.pos1.dP_media_vacuum=u.dP./(x.pos1.alpha*x.pos1.c*x.pos1.V_filt_final/p.A+p.Rm)*p.Rm; 
            % 
            if x.pos1.filtration_finished == 0
                x.pos1.filtration_duration=x.pos1.visc_liq*x.pos1.alpha*x.pos1.c*...
                    (x.pos1.V_filt_final^2-x.pos1.V_filt^2)/(2*p.A^2*u.dP)+...
                    x.pos1.visc_liq*p.Rm*(x.pos1.V_filt_final-x.pos1.V_filt)...
                    /(p.A*u.dP)+x.pos1.filtration_time;     
                residual_filtration = max(x.pos1.filtration_duration - x.pos1.filtration_time,0);
               if residual_filtration < Dt
                    filtration_duration_step=residual_filtration;
                    x.pos1.filtration_finished=1;
                else
                    filtration_duration_step=Dt;
                end
                [x,y]=model_filtration(t,filtration_duration_step,p,u,x,y,n_rotation,1); 
            else
                filtration_duration_step =0;
            end
            deliquoring_duration = max(Dt-filtration_duration_step,0);
            if deliquoring_duration > 0
                x.pos1.rho_liq=p.rho_liquid_phase_from_mass_fr(mean(x.pos1.liq_mass_fr_vect,2));
                x.pos1.S_inf=trapz(x.pos1.x,0.155*(1+0.031*p.N_cap_CSD(x.pos1.x,...
                x.pos1.rho_liq,x.pos1.E,u.dP,x.pos1.L_cake).^(-0.49)).*x.pos1.CSD/x.pos1.m0);              
                [x,y]=model_deliquoring_grad(t+filtration_duration_step,deliquoring_duration,p,u,x,y,n_rotation,1); 
            end
            
            % filtrate volume sensor
            y.cont_sign.pos1_4.t=[y.cont_sign.pos1_4.t y.cont_sign.pos1_4.t(end)-t+...
            y.pos1.(['cycle_' num2str(n_rotation)]).t_filt((end-p.control_time/p.V_filt_step):(end))];
            collected_volume=y.pos1.(['cycle_' num2str(n_rotation)]).V_filt((end-p.control_time/p.V_filt_step):(end));
        end  
        
        %% Position 2
        if n_rotation > 1
            x.pos2.dP_media_vacuum=u.dP./(x.pos2.alpha*x.pos2.c*x.pos2.V_filt_final/p.A+p.Rm)*p.Rm; 
            if x.pos2.filtration_finished == 0
                x.pos2.filtration_duration=x.pos2.visc_liq*x.pos2.alpha*x.pos2.c*...
                    (x.pos2.V_filt_final^2-x.pos2.V_filt^2)/(2*p.A^2*u.dP)+...
                    x.pos2.visc_liq*p.Rm*(x.pos2.V_filt_final-x.pos2.V_filt)...
                    /(p.A*u.dP)+x.pos2.filtration_time;
                residual_filtration = max(x.pos2.filtration_duration - x.pos2.filtration_time,0);
                if residual_filtration < Dt
                    filtration_duration_step=residual_filtration;
                    x.pos2.filtration_finished=1;
                else
                    filtration_duration_step=Dt;
                end
                [x,y]=model_filtration(t,filtration_duration_step,p,u,x,y,n_rotation-1,2); 
            else
                filtration_duration_step=0;
            end
            if x.pos2.washing_finished == 0 && x.pos2.filtration_finished == 1      
                x.pos2.Q_wash=p.A*u.dP/(p.visc_liquid_phase_from_mass_fr(x.pos2.T,[0 1 0]')*(x.pos2.alpha*p.rho_sol*x.pos2.L_cake*(1-x.pos2.E)+p.Rm));          
                x.pos2.washing_duration=(x.pos2.V_wash_final-x.pos2.V_wash)/x.pos2.Q_wash+x.pos2.washing_time;            
                residual_washing = max(x.pos2.washing_duration - x.pos2.washing_time,0); 
                
                if residual_washing < Dt-filtration_duration_step
                    washing_duration_step=residual_washing;
                    x.pos2.washing_finished=1;
                else
                    washing_duration_step=Dt-filtration_duration_step;
                end
                
                [x,y]=model_washing(y.pos2.(['cycle_' num2str(n_rotation-1)]).t_filt(end),...
                    washing_duration_step,p,u,x,y,n_rotation-1,2); 
            else
                washing_duration_step=0;
            end
            deliquoring_duration = max(Dt-filtration_duration_step-washing_duration_step,0);
            if deliquoring_duration > 0
                x.pos2.rho_liq=p.rho_liquid_phase_from_mass_fr(mean(x.pos2.liq_mass_fr_vect,2));
                x.pos2.S_inf=trapz(x.pos2.x,0.155*(1+0.031*p.N_cap_CSD(x.pos2.x,x.pos2.rho_liq,x.pos2.E,u.dP,x.pos2.L_cake).^(-0.49)).*x.pos2.CSD/x.pos2.m0);    % irreducible cake saturation (minimum saturation that can be achieved by displacement of the interstitial liquid by the applied vacuum    
               [x,y]=model_deliquoring_species_grad(t+filtration_duration_step+washing_duration_step,deliquoring_duration,p,u,x,y,n_rotation-1,2); 

            end
            % filtrate volume sensor
            collected_volume=collected_volume+interp1(y.pos2.(['cycle_' num2str(n_rotation-1)]).t_filt,...
                y.pos2.(['cycle_' num2str(n_rotation-1)]).V_filt,...
                y.pos1.(['cycle_' num2str(n_rotation)]).t_filt((end-p.control_time/p.V_filt_step):(end)));
        end
        
        %% Position 3
        if n_rotation > 2
            x.pos3.dP_media_vacuum=u.dP./(x.pos3.alpha*x.pos3.c*x.pos3.V_filt_final/p.A+p.Rm)*p.Rm; 
            if x.pos3.filtration_finished == 0
                x.pos3.filtration_duration=x.pos3.visc_liq*x.pos3.alpha*x.pos3.c*...
                    (x.pos3.V_filt_final^2-x.pos3.V_filt^2)/(2*p.A^2*u.dP)+...
                    x.pos3.visc_liq*p.Rm*(x.pos3.V_filt_final-x.pos3.V_filt)...
                    /(p.A*u.dP)+x.pos3.filtration_time;
                residual_filtration = max(x.pos3.filtration_duration - x.pos3.filtration_time,0);
                if residual_filtration < Dt
                    filtration_duration_step=residual_filtration;
                    x.pos3.filtration_finished=1;
                else
                    filtration_duration_step=Dt;
                end
                [x,y]=model_filtration(t,filtration_duration_step,p,u,x,y,n_rotation-2,3); 
            else
                filtration_duration_step = 0;
            end
            if x.pos3.washing_finished == 0 && x.pos3.filtration_finished == 1 &&...
                    x.pos3.washing_time > 0
                x.pos3.Q_wash=p.A*u.dP/(p.visc_liquid_phase_from_mass_fr(x.pos3.T,[0 1 0]')*(x.pos3.alpha*p.rho_sol*x.pos3.L_cake*(1-x.pos3.E)+p.Rm));          
                x.pos3.washing_duration=(x.pos3.V_wash_final-x.pos3.V_wash)/x.pos3.Q_wash+x.pos3.washing_time;                          
                residual_washing = max(x.pos3.washing_duration - x.pos3.washing_time,0);
                if residual_washing < Dt-filtration_duration_step
                    washing_duration_step=residual_washing;
                    x.pos3.washing_finished=1;
                else
                    washing_duration_step=Dt-filtration_duration_step;
                end
                [x,y]=model_washing(y.pos3.(['cycle_' num2str(n_rotation-2)]).t_filt(end),washing_duration_step,p,u,x,y,n_rotation-2,3); 
            else
                washing_duration_step=0;
            end
            deliquoring_duration = max(Dt-filtration_duration_step-washing_duration_step,0);       
            if deliquoring_duration > 0
                x.pos3.rho_liq=p.rho_liquid_phase_from_mass_fr(mean(x.pos3.liq_mass_fr_vect,2));
                x.pos3.S_inf=trapz(x.pos3.x,0.155*(1+0.031*p.N_cap_CSD(x.pos3.x,x.pos3.rho_liq,x.pos3.E,u.dP,x.pos3.L_cake).^(-0.49)).*x.pos3.CSD/x.pos3.m0);    % irreducible cake saturation (minimum saturation that can be achieved by displacement of the interstitial liquid by the applied vacuum    
                [x,y]=model_deliquoring_species_grad(t+filtration_duration_step+washing_duration_step,deliquoring_duration,p,u,x,y,n_rotation-2,3); 
            end
            % filtrate volume sensor
            collected_volume = collected_volume + interp1(y.pos3.(['cycle_' num2str(n_rotation-2)]).t_filt((end-p.control_time/p.V_filt_step):(end)),...
                y.pos3.(['cycle_' num2str(n_rotation-2)]).V_filt((end-p.control_time/p.V_filt_step):(end)),...
                y.pos1.(['cycle_' num2str(n_rotation)]).t_filt((end-p.control_time/p.V_filt_step):(end)));
        end
        
        
        %% Position 4
        if n_rotation > 3
            x.pos4.dP_media_vacuum=0; % no filter in position 5
            if x.pos4.filtration_finished == 0
                x.pos4.filtration_duration=x.pos4.visc_liq*x.pos4.alpha*x.pos4.c*...
                    (x.pos4.V_filt_final^2-x.pos4.V_filt^2)/(2*p.A^2*u.dP)+...
                    x.pos4.visc_liq*p.Rm*(x.pos4.V_filt_final-x.pos4.V_filt)...
                    /(p.A*u.dP)+x.pos4.filtration_time;
                residual_filtration = max(x.pos4.filtration_duration - x.pos4.filtration_time,0);
                if residual_filtration < Dt
                    filtration_duration_step=residual_filtration;
                    x.pos4.filtration_finished=1;
                else
                    filtration_duration_step=Dt;
                end
                [x,y]=model_filtration(t,filtration_duration_step,p,u,x,y,n_rotation-3,4); 
            else
                filtration_duration_step = 0;
            end
            if x.pos4.washing_finished == 0 && x.pos4.filtration_finished == 1  &&...
                    x.pos4.washing_time > 0
                x.pos2.washing_duration=(x.pos4.V_wash_final-x.pos4.V_wash)/x.pos4.Q_wash+x.pos2.washing_time;                     
                residual_washing = max(x.pos4.washing_duration - x.pos4.washing_time,0);      
                if residual_washing < Dt-filtration_duration_step
                    washing_duration_step=residual_washing;
                    x.pos4.washing_finished=1;
                else
                    washing_duration_step=Dt-filtration_duration_step;
                end
                [x,y]=model_washing(y.pos2.(['cycle_' num2str(n_rotation-3)]).t_filt(end),washing_duration_step,p,u,x,y,n_rotation-3,4);     
            else
                washing_duration_step=0;
            end 
            x.pos4.rho_liq=p.rho_liquid_phase_from_mass_fr(mean(x.pos4.liq_mass_fr_vect,2));
            x.pos4.S_inf=trapz(x.pos4.x,0.155*(1+0.031*p.N_cap_CSD(x.pos4.x,x.pos4.rho_liq,x.pos4.E,u.dP,x.pos4.L_cake).^(-0.49)).*x.pos4.CSD/x.pos4.m0);    % irreducible cake saturation (minimum saturation that can be achieved by displacement of the interstitial liquid by the applied vacuum              
            drying_duration = max(Dt-filtration_duration_step-washing_duration_step,0);       
            if drying_duration > 0
                [x,y]=model_drying(t+filtration_duration_step+washing_duration_step,drying_duration,p,u,x,y,n_rotation-3,4); 
            end
        end    
        
        % filtrate volume sensor
        % collected volume of filtrate in positions 1-4      
        if n_rotation>0
            y.cont_sign.pos1_4.V=[y.cont_sign.pos1_4.V, collected_volume];   
            y.cont_sign.pos1_4.t(isnan(y.cont_sign.pos1_4.V))=[];
            y.cont_sign.pos1_4.V(isnan(y.cont_sign.pos1_4.V))=[];
            
        end
        
        

end
