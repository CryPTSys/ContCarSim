function [x,y]=carousel_simulator(cycle_time,simulation_step,p,d,u,x,y,n_cycle)        
    
        %% Station 1
        if sum(p.ports_working==1)>0
            
            % update pressure drop through filter mesh (depends on u.dP and p.Rm)
            x.pos1.dP_mesh=u.dP./(x.pos1.alpha*x.pos1.c*x.pos1.V_filt_final/p.A+p.Rm(1))*p.Rm(1); 
            
            %----------------------------------------------------------------------------------------
            % filtration 
            if x.pos1.filtration_finished == 0 % filtration not finished
                
               % calculate residual filtration duration
               x.pos1.filtration_duration=x.pos1.visc_liq*x.pos1.alpha*x.pos1.c*... 
                    (x.pos1.V_filt_final^2-x.pos1.V_filt^2)/(2*p.A^2*u.dP)+...
                    x.pos1.visc_liq*p.Rm(1)*(x.pos1.V_filt_final-x.pos1.V_filt)...
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
        else
            y.pos1.(['cycle_' num2str(n_cycle)]).t_filt=[y.pos1.(['cycle_' num2str(n_cycle)]).t_filt cycle_time:p.filtration_sampling_time:cycle_time+simulation_step];
            y.pos1.(['cycle_' num2str(n_cycle)]).V_filt=[y.pos1.(['cycle_' num2str(n_cycle)]).V_filt zeros(1, length(cycle_time:p.filtration_sampling_time:cycle_time+simulation_step))];
            
        end
        
        % filtrate volume sensor WI101 - add filtrate volume from Station 1
        sampling_times1_4=y.pos1.(['cycle_' num2str(n_cycle)]).t_filt((end+1-p.control_interval/p.filtration_sampling_time):end);
        y.cont_sign.pos1_4.t=[y.cont_sign.pos1_4.t y.cont_sign.pos1_4.t(end)-cycle_time+sampling_times1_4];
        collected_volume=y.pos1.(['cycle_' num2str(n_cycle)]).V_filt((end+1-p.control_interval/p.filtration_sampling_time):end)+y.cont_sign.pos1_4.bias;
        
        %% Station 2
        if sum(p.ports_working==2)>0
            % update pressure drop through filter mesh (depends on u.dP and p.Rm)
            x.pos2.dP_mesh=u.dP./(x.pos2.alpha*x.pos2.c*x.pos2.V_filt_final/p.A+p.Rm(2))*p.Rm(2); 
            
            %----------------------------------------------------------------------------------------
            % filtration
            
            if x.pos2.filtration_finished == 0 % filtration not finished
                
                % calculate residual filtration duration 
                x.pos2.filtration_duration=x.pos2.visc_liq*x.pos2.alpha*x.pos2.c*...
                    (x.pos2.V_filt_final^2-x.pos2.V_filt^2)/(2*p.A^2*u.dP)+...
                    x.pos2.visc_liq*p.Rm(2)*(x.pos2.V_filt_final-x.pos2.V_filt)...
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
            % deliquoring            
            
            % calculate duration of deliquoring step in current simulation step
            deliquoring_duration = max(simulation_step-filtration_duration_step,0);
            
            % simulate deliquoring
            if deliquoring_duration > 0
               [x,y]=model_deliquoring_grad(cycle_time+filtration_duration_step,deliquoring_duration,p,u,x,y,n_cycle-1,2); 
            end
            %----------------------------------------------------------------------------------------
            % filtrate volume sensor WI101 - add filtrate volume from Station 2
            collected_volume=collected_volume+interp1(y.pos2.(['cycle_' num2str(n_cycle-1)]).t_filt,...
                y.pos2.(['cycle_' num2str(n_cycle-1)]).V_filt,...
                sampling_times1_4);
            
        elseif n_cycle>1
            y.pos2.(['cycle_' num2str(n_cycle-1)]).t_filt=[y.pos2.(['cycle_' num2str(n_cycle-1)]).t_filt cycle_time:p.filtration_sampling_time:cycle_time+simulation_step];
            y.pos2.(['cycle_' num2str(n_cycle-1)]).V_filt=[y.pos2.(['cycle_' num2str(n_cycle-1)]).V_filt zeros(1, length(cycle_time:p.filtration_sampling_time:cycle_time+simulation_step))];
            

        end
        
        %% Station 3
        if sum(p.ports_working==3)>0
            % update pressure drop through filter mesh (depends on u.dP and p.Rm)
            x.pos3.dP_mesh=u.dP./(x.pos3.alpha*x.pos3.c*x.pos3.V_filt_final/p.A+p.Rm(3))*p.Rm(3); 
            
            %----------------------------------------------------------------------------------------
            % filtration
            
            if x.pos3.filtration_finished == 0 % filtration not finished
                
                % calculate residual filtration duration 
                x.pos3.filtration_duration=x.pos3.visc_liq*x.pos3.alpha*x.pos3.c*...
                    (x.pos3.V_filt_final^2-x.pos3.V_filt^2)/(2*p.A^2*u.dP)+...
                    x.pos3.visc_liq*p.Rm(3)*(x.pos3.V_filt_final-x.pos3.V_filt)...
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
            % deliquoring            
            
            % calculate duration of deliquoring step in current simulation step
            deliquoring_duration = max(simulation_step-filtration_duration_step,0);   
            
            % simulate deliquoring
            if deliquoring_duration > 0
               [x,y]=model_deliquoring_grad(cycle_time+filtration_duration_step,deliquoring_duration,p,u,x,y,n_cycle-2,3); 
            end
                        
            %----------------------------------------------------------------------------------------
            % filtrate volume sensor WI101 - add filtrate volume from Station 3
            collected_volume = collected_volume + interp1(y.pos3.(['cycle_' num2str(n_cycle-2)]).t_filt,...
                y.pos3.(['cycle_' num2str(n_cycle-2)]).V_filt,...
                sampling_times1_4);
            
        elseif n_cycle > 2
            y.pos3.(['cycle_' num2str(n_cycle-2)]).t_filt=[y.pos3.(['cycle_' num2str(n_cycle-2)]).t_filt cycle_time:p.filtration_sampling_time:cycle_time+simulation_step];
            y.pos3.(['cycle_' num2str(n_cycle-2)]).V_filt=[y.pos3.(['cycle_' num2str(n_cycle-2)]).V_filt zeros(1, length(cycle_time:p.filtration_sampling_time:cycle_time+simulation_step))];
            
        end        
        
        %% Station 4
        if sum(p.ports_working==4)>0
            % update pressure drop through filter mesh (depends on u.dP_drying and p.Rm)
            x.pos4.dP_mesh=u.dP_drying./(x.pos4.alpha*x.pos4.c*x.pos4.V_filt_final/p.A+p.Rm(4))*p.Rm(4); 
            
            %----------------------------------------------------------------------------------------
            % filtration            
            if x.pos4.filtration_finished == 0 % filtration not finished
                
                % calculate residual filtration duration 
                x.pos4.filtration_duration=x.pos4.visc_liq*x.pos4.alpha*x.pos4.c*...
                    (x.pos4.V_filt_final^2-x.pos4.V_filt^2)/(2*p.A^2*u.dP)+...
                    x.pos4.visc_liq*p.Rm(4)*(x.pos4.V_filt_final-x.pos4.V_filt)...
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
            % drying
            
            % calculate duration of drying step in current simulation step
            drying_duration = max(simulation_step-filtration_duration_step,0);   
%             x.pos4.S=x.pos4.S./x.pos4.S*0.10;
            % drying simulation
            if drying_duration > 0
                [x,y]=model_drying(cycle_time+filtration_duration_step,drying_duration,p,d,u,x,y,n_cycle-3,4); 
            end
            
            %----------------------------------------------------------------------------------------
            % filtrate volume sensor WI101 - add filtrate volume from Station 4 (if there is any) 
            if sum(y.pos4.(['cycle_' num2str(n_cycle-3)]).V_filt)>0
                collected_volume = collected_volume + interp1(y.pos4.(['cycle_' num2str(n_cycle-3)]).t_filt,...
                    y.pos4.(['cycle_' num2str(n_cycle-3)]).V_filt,...
                    sampling_times1_4);
            end
        elseif n_cycle > 3
            y.pos4.(['cycle_' num2str(n_cycle-3)]).t_filt=[y.pos4.(['cycle_' num2str(n_cycle-3)]).t_filt cycle_time:p.filtration_sampling_time:cycle_time+simulation_step];
            y.pos4.(['cycle_' num2str(n_cycle-3)]).V_filt=[y.pos4.(['cycle_' num2str(n_cycle-3)]).V_filt zeros(1, length(cycle_time:p.filtration_sampling_time:cycle_time+simulation_step))];            
            y.pos4.(['cycle_' num2str(n_cycle-3)]).t_drying=[y.pos4.(['cycle_' num2str(n_cycle-3)]).t_drying cycle_time:p.filtration_sampling_time:cycle_time+simulation_step];
            y.pos4.(['cycle_' num2str(n_cycle-3)]).Tg=[y.pos4.(['cycle_' num2str(n_cycle-3)]).Tg ones(1, length(cycle_time:p.filtration_sampling_time:cycle_time+simulation_step))*295.15];
        end      
        
        %------------------------------------------------------------------
        % update filtrate sensor WI101 reading
        y.cont_sign.pos1_4.V=[y.cont_sign.pos1_4.V, collected_volume];   
        y.cont_sign.pos1_4.t(isnan(y.cont_sign.pos1_4.V))=[];
        y.cont_sign.pos1_4.V(isnan(y.cont_sign.pos1_4.V))=[];            
               
end