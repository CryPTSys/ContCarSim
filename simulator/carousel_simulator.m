function [x,y,measurements, measurements_nf]=carousel_simulator(cycle_time,simulation_step,p,d,u,x,y,measurements,measurements_nf,n_cycle)        
        %% Station 1
        if p.stations_working(1)==1

            % update pressure drop through filter mesh (depends on u.P_compr and p.Rm)
            x.pos1.P_compr_mesh=u.P_compr./(x.pos1.alpha*x.pos1.c*x.pos1.V_filt_final/p.A+p.Rm(1))*p.Rm(1); 
            
            %----------------------------------------------------------------------------------------
            % filtration 
            if x.pos1.filtration_finished == 0 % filtration not finished
            
               % calculate residual filtration duration
               x.pos1.filtration_duration=x.pos1.visc_liq*x.pos1.alpha*x.pos1.c*... 
                    (x.pos1.V_filt_final^2-x.pos1.V_filt^2)/(2*p.A^2*u.P_compr)+...
                    x.pos1.visc_liq*p.Rm(1)*(x.pos1.V_filt_final-x.pos1.V_filt)...
                    /(p.A*u.P_compr)+x.pos1.filtration_time;     
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
            y.pos1.(['batch_' num2str(n_cycle)]).t=[y.pos1.(['batch_' num2str(n_cycle)]).t cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step];
            y.pos1.(['batch_' num2str(n_cycle)]).m_filt=[y.pos1.(['batch_' num2str(n_cycle)]).m_filt zeros(1, length(cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step))];
            
        end
                
            
        
        % filtrate volume sensor WI101 - add filtrate volume from Station 1
        sampling_intervals1_4=y.pos1.(['batch_' num2str(n_cycle)]).t((end+1-p.integration_interval/p.filtration_sampling_interval):end);
        measurements.t_meas=[measurements.t_meas measurements.t_meas(end)-cycle_time+sampling_intervals1_4];
        collected_mass=y.pos1.(['batch_' num2str(n_cycle)]).m_filt((end+1-p.integration_interval/p.filtration_sampling_interval):end)+x.m_filt_bias;
        
        measurements_nf.t_meas=[measurements_nf.t_meas measurements_nf.t_meas(end)-cycle_time+sampling_intervals1_4];
        %% Station 2
        if p.stations_working(2)==1
            % update pressure drop through filter mesh (depends on u.P_compr and p.Rm)
            x.pos2.P_compr_mesh=u.P_compr./(x.pos2.alpha*x.pos2.c*x.pos2.V_filt_final/p.A+p.Rm(2))*p.Rm(2); 
            
            %----------------------------------------------------------------------------------------
            % filtration
            
            if x.pos2.filtration_finished == 0 % filtration not finished
                
                % calculate residual filtration duration 
                x.pos2.filtration_duration=x.pos2.visc_liq*x.pos2.alpha*x.pos2.c*...
                    (x.pos2.V_filt_final^2-x.pos2.V_filt^2)/(2*p.A^2*u.P_compr)+...
                    x.pos2.visc_liq*p.Rm(2)*(x.pos2.V_filt_final-x.pos2.V_filt)...
                    /(p.A*u.P_compr)+x.pos2.filtration_time;
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
            collected_mass=collected_mass+interp1(y.pos2.(['batch_' num2str(n_cycle-1)]).t,...
                y.pos2.(['batch_' num2str(n_cycle-1)]).m_filt,...
                sampling_intervals1_4);
            
        elseif n_cycle>1
            y.pos2.(['batch_' num2str(n_cycle-1)]).t=[y.pos2.(['batch_' num2str(n_cycle-1)]).t cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step];
            y.pos2.(['batch_' num2str(n_cycle-1)]).m_filt=[y.pos2.(['batch_' num2str(n_cycle-1)]).m_filt zeros(1, length(cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step))];
            

        end
        
        %% Station 3
        if p.stations_working(3)==1
            % update pressure drop through filter mesh (depends on u.P_compr and p.Rm)
            x.pos3.P_compr_mesh=u.P_compr./(x.pos3.alpha*x.pos3.c*x.pos3.V_filt_final/p.A+p.Rm(3))*p.Rm(3); 
            
            %----------------------------------------------------------------------------------------
            % filtration
            
            if x.pos3.filtration_finished == 0 % filtration not finished
                
                % calculate residual filtration duration 
                x.pos3.filtration_duration=x.pos3.visc_liq*x.pos3.alpha*x.pos3.c*...
                    (x.pos3.V_filt_final^2-x.pos3.V_filt^2)/(2*p.A^2*u.P_compr)+...
                    x.pos3.visc_liq*p.Rm(3)*(x.pos3.V_filt_final-x.pos3.V_filt)...
                    /(p.A*u.P_compr)+x.pos3.filtration_time;
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
            collected_mass = collected_mass + interp1(y.pos3.(['batch_' num2str(n_cycle-2)]).t,...
                y.pos3.(['batch_' num2str(n_cycle-2)]).m_filt,...
                sampling_intervals1_4);
            
        elseif n_cycle > 2
            y.pos3.(['batch_' num2str(n_cycle-2)]).t=[y.pos3.(['batch_' num2str(n_cycle-2)]).t cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step];
            y.pos3.(['batch_' num2str(n_cycle-2)]).m_filt=[y.pos3.(['batch_' num2str(n_cycle-2)]).m_filt zeros(1, length(cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step))];
            
        end        
        
        %% Station 4
        if p.stations_working(4)==1
            % update pressure drop through filter mesh (depends on u.P_compr_drying and p.Rm)
            x.pos4.P_compr_mesh=u.P_compr./(x.pos4.alpha*x.pos4.c*x.pos4.V_filt_final/p.A+p.Rm(4))*p.Rm(4); 
            
            %----------------------------------------------------------------------------------------
            % filtration            
            if x.pos4.filtration_finished == 0 % filtration not finished
                
                % calculate residual filtration duration 
                x.pos4.filtration_duration=x.pos4.visc_liq*x.pos4.alpha*x.pos4.c*...
                    (x.pos4.V_filt_final^2-x.pos4.V_filt^2)/(2*p.A^2*u.P_compr)+...
                    x.pos4.visc_liq*p.Rm(4)*(x.pos4.V_filt_final-x.pos4.V_filt)...
                    /(p.A*u.P_compr)+x.pos4.filtration_time;
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
            else
                y.pos4.(['batch_' num2str(n_cycle-3)]).Tg_top=[y.pos4.(['batch_' num2str(n_cycle-3)]).Tg_top ones(1, length(cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step))*p.T_room];
                y.pos4.(['batch_' num2str(n_cycle-3)]).Tg_bot=[y.pos4.(['batch_' num2str(n_cycle-3)]).Tg_bot ones(1, length(cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step))*p.T_room];
                y.pos4.(['batch_' num2str(n_cycle-3)]).Ts_bot=[y.pos4.(['batch_' num2str(n_cycle-3)]).Ts_bot ones(1, length(cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step))*p.T_room];
                y.pos4.(['batch_' num2str(n_cycle-3)]).w_EtOH_gas=[y.pos4.(['batch_' num2str(n_cycle-3)]).w_EtOH_gas zeros(1, length(cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step))];
                y.pos4.(['batch_' num2str(n_cycle-3)]).Vdryer=[y.pos4.(['batch_' num2str(n_cycle-3)]).Vdryer zeros(1, length(cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step))];
            end
            
            %----------------------------------------------------------------------------------------
            % filtrate volume sensor WI101 - add filtrate volume from Station 4 (if there is any) 
            if sum(y.pos4.(['batch_' num2str(n_cycle-3)]).m_filt)>0
                collected_mass = collected_mass + interp1(y.pos4.(['batch_' num2str(n_cycle-3)]).t(1:length(y.pos4.(['batch_' num2str(n_cycle-3)]).m_filt)),...
                    y.pos4.(['batch_' num2str(n_cycle-3)]).m_filt,...
                    sampling_intervals1_4);
            end
        elseif n_cycle > 3
            y.pos4.(['batch_' num2str(n_cycle-3)]).t=[y.pos4.(['batch_' num2str(n_cycle-3)]).t cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step];
            y.pos4.(['batch_' num2str(n_cycle-3)]).m_filt=[y.pos4.(['batch_' num2str(n_cycle-3)]).m_filt zeros(1, length(cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step))];                        
            y.pos4.(['batch_' num2str(n_cycle-3)]).Tg_top=[y.pos4.(['batch_' num2str(n_cycle-3)]).Tg_top ones(1, length(cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step))*p.T_room];
            y.pos4.(['batch_' num2str(n_cycle-3)]).Tg_bot=[y.pos4.(['batch_' num2str(n_cycle-3)]).Tg_bot ones(1, length(cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step))*p.T_room];
            y.pos4.(['batch_' num2str(n_cycle-3)]).Ts_bot=[y.pos4.(['batch_' num2str(n_cycle-3)]).Ts_bot ones(1, length(cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step))*p.T_room];
            y.pos4.(['batch_' num2str(n_cycle-3)]).w_EtOH_gas=[y.pos4.(['batch_' num2str(n_cycle-3)]).w_EtOH_gas zeros(1, length(cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step))];
            y.pos4.(['batch_' num2str(n_cycle-3)]).w_EtOH_cake=[y.pos4.(['batch_' num2str(n_cycle-3)]).w_EtOH_cake zeros(1, length(cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step))];
            y.pos4.(['batch_' num2str(n_cycle-3)]).Vdryer=[y.pos4.(['batch_' num2str(n_cycle-3)]).Vdryer zeros(1, length(cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step))];
            y.pos4.(['batch_' num2str(n_cycle-3)]).S=[y.pos4.(['batch_' num2str(n_cycle-3)]).S zeros(1, length(cycle_time:p.filtration_sampling_interval:cycle_time+simulation_step))];
        end      

        %------------------------------------------------------------------
        % Sensors
        %------------------------------------------------------------------
        % filtrate sensor WI101 reading
        measurements.m_filt_WI101=[measurements.m_filt_WI101, collected_mass+randn(1,length(collected_mass))*0.00005];   
        measurements.t_meas(isnan(measurements.m_filt_WI101))=[];
        measurements.m_filt_WI101(isnan(measurements.m_filt_WI101))=[];  
        
        measurements_nf.m_filt_WI101=[measurements_nf.m_filt_WI101, collected_mass];   
        measurements_nf.t_meas(isnan(measurements_nf.m_filt_WI101))=[];
        measurements_nf.m_filt_WI101(isnan(measurements_nf.m_filt_WI101))=[];  
        % slurry concentration reading AI101
        measurements.c_slurry_AI101=[measurements.c_slurry_AI101 x.pos1.c_slurry*ones(1,length(collected_mass))+randn(1,length(collected_mass))*0.05];
        measurements.c_slurry_AI101(measurements.c_slurry_AI101<1)=0;
        
        measurements_nf.c_slurry_AI101=[measurements_nf.c_slurry_AI101, ...
            x.pos1.c_slurry*ones(1,length(collected_mass))];
        % camera vision system LI101
        measurements.L_cake_LI101=[measurements.L_cake_LI101 x.pos1.L_cake*ones(1,length(collected_mass))];
        measurements_nf.L_cake_LI101=[measurements_nf.L_cake_LI101 x.pos1.L_cake*ones(1,length(collected_mass))];
        
        measurements.V_slurry_LI101=[measurements.V_slurry_LI101 x.pos1.V_slurry_initial*ones(1,length(collected_mass))];
        measurements_nf.V_slurry_LI101=[measurements_nf.V_slurry_LI101 x.pos1.V_slurry_initial*ones(1,length(collected_mass))];
        
        % pressure sensor PI101
        measurements.P_PI101=[measurements.P_PI101 u.P_compr*ones(1,length(collected_mass))];
        
        measurements_nf.P_PI101=[measurements_nf.P_PI101 u.P_compr*ones(1,length(collected_mass))]; 
        
        % pressure sensor PI102
        measurements.P_PI102=[measurements.P_PI102 0*ones(1,length(collected_mass))];
        
        measurements_nf.P_PI102=[measurements_nf.P_PI102 0*ones(1,length(collected_mass))];
        
        % temperature sensor TI101
        measurements.Tg_in_TI101=[measurements.Tg_in_TI101 round(u.Tinlet_drying,1)*ones(1,length(collected_mass))];
        
        measurements_nf.Tg_in_TI101=[measurements_nf.Tg_in_TI101 u.Tinlet_drying*ones(1,length(collected_mass))];
        
        % temperature sensor TI102
        if n_cycle <= 3
            measurements.Tg_out_TI102=[measurements.Tg_out_TI102, ...
                measurements.Tg_out_TI102(end)*ones(1,length(collected_mass))];
            
            measurements_nf.Tg_out_TI102=[measurements_nf.Tg_out_TI102, ...
                measurements_nf.Tg_out_TI102(end)*ones(1,length(collected_mass))];
        else
            measurements.Tg_out_TI102=[measurements.Tg_out_TI102, ...
                round(y.pos4.(['batch_' num2str(n_cycle-3)]).Tg_bot((end+1-...
                p.integration_interval/p.filtration_sampling_interval):end),1)];
            
            measurements_nf.Tg_out_TI102=[measurements_nf.Tg_out_TI102, ...
                y.pos4.(['batch_' num2str(n_cycle-3)]).Tg_bot((end+1-...
                p.integration_interval/p.filtration_sampling_interval):end)];
        end
        % flowrate sensor FI101
        if n_cycle <= 3
            measurements.Vdryer_FI101=[measurements.Vdryer_FI101, ...
                zeros(1,length(collected_mass))];
            measurements_nf.Vdryer_FI101=[measurements_nf.Vdryer_FI101, ...
                zeros(1,length(collected_mass))];
        else
            measurements.Vdryer_FI101=[measurements.Vdryer_FI101, ...
                round(y.pos4.(['batch_' num2str(n_cycle-3)]).Vdryer((end+1-...
                p.integration_interval/p.filtration_sampling_interval):end),1)];
            measurements_nf.Vdryer_FI101=[measurements_nf.Vdryer_FI101, ...
                y.pos4.(['batch_' num2str(n_cycle-3)]).Vdryer((end+1-...
                p.integration_interval/p.filtration_sampling_interval):end)];
        end
end