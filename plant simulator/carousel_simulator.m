function [x,y]=carousel_simulator(cycle_time,simulation_step,p,u,x,y,n_cycle)
        
        %% Slurry tank
%         x.slurry_tank_volume = x.slurry_tank_volume + ...
%             cryst_output.flowrate_MSMPR*simulation_step - u.flowrate_slurry*simulation_step;
    
        %% Position 0 - charge cell
        x.pos0.charge_cell_volume = u.V_slurry;%;x.pos0.charge_cell_volume + u.V_slurry;
        
        %% Position 1: filtration (+ pre-deliquoring_grad)
        if n_cycle > 0
            
            % update pressure drop through filter mesh (depends on u.dP and p.Rm)
            x.pos1.dP_media_vacuum=u.dP./(x.pos1.alpha*x.pos1.c*x.pos1.V_filt_final/p.A+p.Rm)*p.Rm; 
            
            %----------------------------------------------------------------------------------------
            % filtration 
            if x.pos1.filtration_finished == 0 % filtration not finished
                
                % calculate residual filtration duration
                x.pos1.filtration_duration=x.pos1.visc_liq*x.pos1.alpha*x.pos1.c*... 
                    (x.pos1.V_filt_final^2-x.pos1.V_filt^2)/(2*p.A^2*u.dP)+...
                    x.pos1.visc_liq*p.Rm*(x.pos1.V_filt_final-x.pos1.V_filt)...
                    /(p.A*u.dP)+x.pos1.filtration_time;     
                residual_filtration = max(x.pos1.filtration_duration - x.pos1.filtration_time,0);
                
               % calculate duration of filtration step in current simulation step
               if residual_filtration < simulation_step
                    filtration_duration_step=residual_filtration;
               else
                    filtration_duration_step=simulation_step;
               end
               
               % filtration simulation
               [x,y]=model_filtration(cycle_time,filtration_duration_step,p,u,x,y,n_cycle,1); 
               
            else % filtration is finished - no filtration simulation
                filtration_duration_step =0;
            end
            
            %----------------------------------------------------------------------------------------
            % deliquoring
            
            % calculate duration of deliquoring step in current simulation step
            deliquoring_duration = max(simulation_step-filtration_duration_step,0);
            
            % simulate deliquoring
            if deliquoring_duration > 0                       
                [x,y]=model_deliquoring_grad(cycle_time+filtration_duration_step,deliquoring_duration,p,u,x,y,n_cycle,1); 
            end
            
            %----------------------------------------------------------------------------------------
            % filtrate volume sensor - add filtrate volume from Station 1
            sampling_times1_4=y.pos1.(['cycle_' num2str(n_cycle)]).t_filt((end-p.control_interval/p.filtration_sampling_time):end);
            y.cont_sign.pos1_4.t=[y.cont_sign.pos1_4.t y.cont_sign.pos1_4.t(end)-cycle_time+sampling_times1_4];
            collected_volume=y.pos1.(['cycle_' num2str(n_cycle)]).V_filt((end-p.control_interval/p.filtration_sampling_time):end);
        end  
        
        %% Position 2
        if n_cycle > 1
            % update pressure drop through filter mesh (depends on u.dP and p.Rm)
            x.pos2.dP_media_vacuum=u.dP./(x.pos2.alpha*x.pos2.c*x.pos2.V_filt_final/p.A+p.Rm)*p.Rm; 
            
            %----------------------------------------------------------------------------------------
            % filtration
            
            if x.pos2.filtration_finished == 0 % filtration not finished
                
                % calculate residual filtration duration 
                x.pos2.filtration_duration=x.pos2.visc_liq*x.pos2.alpha*x.pos2.c*...
                    (x.pos2.V_filt_final^2-x.pos2.V_filt^2)/(2*p.A^2*u.dP)+...
                    x.pos2.visc_liq*p.Rm*(x.pos2.V_filt_final-x.pos2.V_filt)...
                    /(p.A*u.dP)+x.pos2.filtration_time;
                residual_filtration = max(x.pos2.filtration_duration - x.pos2.filtration_time,0);
                
                % calculate duration of filtration step in current simulation step
                if residual_filtration < simulation_step
                    filtration_duration_step=residual_filtration;
                else
                    filtration_duration_step=simulation_step;
                end
                
                % filtration simulation
                [x,y]=model_filtration(cycle_time,filtration_duration_step,p,u,x,y,n_cycle-1,2); 
                
            else % filtration is finished - no filtration simulation
                filtration_duration_step=0;
            end
            
            
            %----------------------------------------------------------------------------------------
            % washing
            
            if x.pos2.washing_finished == 0 && x.pos2.filtration_finished == 1 % washing not finished, filtration finished     
                
                % calculate residual washing duration
                x.pos2.Q_wash=p.A*u.dP/(p.visc_liquid_phase_from_mass_fr(x.pos2.T,p.wash_solvent_mass_fr)*(x.pos2.alpha*p.rho_sol*x.pos2.L_cake*(1-x.pos2.E)+p.Rm)); % wash solvent flowrate with current dP         
                x.pos2.washing_duration=(x.pos2.V_wash_final-x.pos2.V_wash)/x.pos2.Q_wash+x.pos2.washing_time; % sums duration of previous washing steps with time needed for completing washing with current u vector            
                residual_washing = max(x.pos2.washing_duration - x.pos2.washing_time,0); % residual washing time
                
                % calculate duration of washing step in current simulation step
                if residual_washing < simulation_step-filtration_duration_step
                    washing_duration_step=residual_washing;
                    x.pos2.washing_finished=1;
                else
                    washing_duration_step=simulation_step-filtration_duration_step;
                end
                
                % washing simulation
                [x,y]=model_washing(y.pos2.(['cycle_' num2str(n_cycle-1)]).t_filt(end),...
                    washing_duration_step,p,u,x,y,n_cycle-1,2); 
            
            else % washing doesn't happen
                washing_duration_step=0;
            end
            
            %----------------------------------------------------------------------------------------
            % deliquoring            
            
            % calculate duration of deliquoring step in current simulation step
            deliquoring_duration = max(simulation_step-filtration_duration_step-washing_duration_step,0);
            
            % simulate deliquoring
            if deliquoring_duration > 0
               [x,y]=model_deliquoring_species_grad(cycle_time+filtration_duration_step+washing_duration_step,deliquoring_duration,p,u,x,y,n_cycle-1,2); 
            end
            
            %----------------------------------------------------------------------------------------
            % filtrate volume sensor - add filtrate volume from Station 2
            collected_volume=collected_volume+interp1(y.pos2.(['cycle_' num2str(n_cycle-1)]).t_filt,...
                y.pos2.(['cycle_' num2str(n_cycle-1)]).V_filt,...
                sampling_times1_4);
        end
        
        %% Position 3
        if n_cycle > 2
            % update pressure drop through filter mesh (depends on u.dP and p.Rm)
            x.pos3.dP_media_vacuum=u.dP./(x.pos3.alpha*x.pos3.c*x.pos3.V_filt_final/p.A+p.Rm)*p.Rm; 
            
            %----------------------------------------------------------------------------------------
            % filtration
            
            if x.pos3.filtration_finished == 0 % filtration not finished
                
                % calculate residual filtration duration 
                x.pos3.filtration_duration=x.pos3.visc_liq*x.pos3.alpha*x.pos3.c*...
                    (x.pos3.V_filt_final^2-x.pos3.V_filt^2)/(2*p.A^2*u.dP)+...
                    x.pos3.visc_liq*p.Rm*(x.pos3.V_filt_final-x.pos3.V_filt)...
                    /(p.A*u.dP)+x.pos3.filtration_time;
                residual_filtration = max(x.pos3.filtration_duration - x.pos3.filtration_time,0);
                
                % calculate duration of filtration step in current simulation step
                if residual_filtration < simulation_step
                    filtration_duration_step=residual_filtration;                    
                else
                    filtration_duration_step=simulation_step;
                end
                
                % filtration simulation
                [x,y]=model_filtration(cycle_time,filtration_duration_step,p,u,x,y,n_cycle-2,3); 
                
            else % filtration is finished - no filtration simulation               
                filtration_duration_step = 0;                
            end
            
            %----------------------------------------------------------------------------------------
            % washing
            
            if x.pos3.washing_finished == 0 && x.pos3.filtration_finished == 1 &&...
                    x.pos3.washing_time > 0 % washing started but not finished, filtration finished  
                
                % calculate residual washing duration
                x.pos3.Q_wash=p.A*u.dP/(p.visc_liquid_phase_from_mass_fr(x.pos3.T,[0 1 0]')*(x.pos3.alpha*p.rho_sol*x.pos3.L_cake*(1-x.pos3.E)+p.Rm));          
                x.pos3.washing_duration=(x.pos3.V_wash_final-x.pos3.V_wash)/x.pos3.Q_wash+x.pos3.washing_time;                          
                residual_washing = max(x.pos3.washing_duration - x.pos3.washing_time,0);
                
                % calculate duration of washing step in current simulation step
                if residual_washing < simulation_step-filtration_duration_step
                    washing_duration_step=residual_washing;
                    x.pos3.washing_finished=1;
                else
                    washing_duration_step=simulation_step-filtration_duration_step;
                end
                
                % washing simulation
                [x,y]=model_washing(y.pos3.(['cycle_' num2str(n_cycle-2)]).t_filt(end),washing_duration_step,p,u,x,y,n_cycle-2,3); 
            
            else % washing doesn't happen
                washing_duration_step=0;
            end
            
            %----------------------------------------------------------------------------------------
            % deliquoring            
            
            % calculate duration of deliquoring step in current simulation step
            deliquoring_duration = max(simulation_step-filtration_duration_step-washing_duration_step,0);   
            
            % simulate deliquoring
            if deliquoring_duration > 0
               [x,y]=model_deliquoring_species_grad(cycle_time+filtration_duration_step+washing_duration_step,deliquoring_duration,p,u,x,y,n_cycle-2,3); 
            end
                        
            %----------------------------------------------------------------------------------------
            % filtrate volume sensor - add filtrate volume from Station 3
            collected_volume = collected_volume + interp1(y.pos3.(['cycle_' num2str(n_cycle-2)]).t_filt((end-p.control_interval/p.filtration_sampling_time):(end)),...
                y.pos3.(['cycle_' num2str(n_cycle-2)]).V_filt((end-p.control_interval/p.filtration_sampling_time):(end)),...
                sampling_times1_4);
        end
        
        
        %% Position 4
        if n_cycle > 3
            % update pressure drop through filter mesh (depends on u.dP_drying and p.Rm)
            x.pos4.dP_media_vacuum=u.dP_drying./(x.pos4.alpha*x.pos4.c*x.pos4.V_filt_final/p.A+p.Rm)*p.Rm; 
            
            %----------------------------------------------------------------------------------------
            % filtration
            
            if x.pos4.filtration_finished == 0 % filtration not finished
                
                % calculate residual filtration duration 
                x.pos4.filtration_duration=x.pos4.visc_liq*x.pos4.alpha*x.pos4.c*...
                    (x.pos4.V_filt_final^2-x.pos4.V_filt^2)/(2*p.A^2*u.dP)+...
                    x.pos4.visc_liq*p.Rm*(x.pos4.V_filt_final-x.pos4.V_filt)...
                    /(p.A*u.dP)+x.pos4.filtration_time;
                residual_filtration = max(x.pos4.filtration_duration - x.pos4.filtration_time,0);
                
                % calculate duration of filtration step in current simulation step
                if residual_filtration < simulation_step
                    filtration_duration_step=residual_filtration;
                else
                    filtration_duration_step=simulation_step;
                end
                
                % filtration simulation
                [x,y]=model_filtration(cycle_time,filtration_duration_step,p,u,x,y,n_cycle-3,4); 
                
            else % filtration is finished - no filtration simulation 
                filtration_duration_step = 0;
            end
            
            %----------------------------------------------------------------------------------------
            % washing
                     
            if x.pos4.washing_finished == 0 && x.pos4.filtration_finished == 1  &&...
                    x.pos4.washing_time > 0 % washing started but not finished, filtration finished  
                
                % calculate residual washing duration
                x.pos4.washing_duration=(x.pos4.V_wash_final-x.pos4.V_wash)/x.pos4.Q_wash+x.pos4.washing_time;                     
                residual_washing = max(x.pos4.washing_duration - x.pos4.washing_time,0); 
                
                % calculate duration of washing step in current simulation step
                if residual_washing < simulation_step-filtration_duration_step
                    washing_duration_step=residual_washing;
                    x.pos4.washing_finished=1;
                else
                    washing_duration_step=simulation_step-filtration_duration_step;
                end
                
                % washing simulation
                [x,y]=model_washing(y.pos4.(['cycle_' num2str(n_cycle-3)]).t_filt(end),washing_duration_step,p,u,x,y,n_cycle-3,4);     
            
            else % washing doesn't happen
                washing_duration_step=0;
            end
            
            %----------------------------------------------------------------------------------------
            % drying
            
            % calculate duration of drying step in current simulation step
            drying_duration = max(simulation_step-filtration_duration_step-washing_duration_step,0);   
            
            % drying simulation
            if drying_duration > 0
                [x,y]=model_drying(cycle_time+filtration_duration_step+washing_duration_step,drying_duration,p,u,x,y,n_cycle-3,4); 
            end
            
            %----------------------------------------------------------------------------------------
            % filtrate volume sensor - add filtrate volume from Station 4 (if there is any) 
            if sum(y.pos4.(['cycle_' num2str(n_cycle-3)]).V_filt)>0
                collected_volume = collected_volume + interp1(y.pos4.(['cycle_' num2str(n_cycle-3)]).t_filt((end-p.control_interval/p.filtration_sampling_time):(end)),...
                    y.pos4.(['cycle_' num2str(n_cycle-3)]).V_filt((end-p.control_interval/p.filtration_sampling_time):(end)),...
                    sampling_times1_4);
            end
        end    
        
        %------------------------------------------------------------------
        % update filtrate volume sensor reading
        if n_cycle>0
            y.cont_sign.pos1_4.V=[y.cont_sign.pos1_4.V, collected_volume];   
            y.cont_sign.pos1_4.t(isnan(y.cont_sign.pos1_4.V))=[];
            y.cont_sign.pos1_4.V(isnan(y.cont_sign.pos1_4.V))=[];            
        end
               
end