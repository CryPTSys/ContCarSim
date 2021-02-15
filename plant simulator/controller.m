function u = controller(sp,t,u_ss,p,u,y,n_rotation, control_flag)
    % split range control - MVs=Dp_drying/Tg,out_drying, CV=Tg,out_drying
    if control_flag > 0
        
        %% error calculation
        if n_rotation>4 % starts from second drying cycle
            
            % re-sampling of set-point temperature profile based on current
            % cycle duration
            control_times_ref=(p.control_interval*sp.t_rot_ref/u.t_rot):...
                (p.control_interval*sp.t_rot_ref/u.t_rot):(sp.t_rot_ref);
            sp.Tg_pos4=interp1(sp.t_ref_Tg,sp.Tg_pos4_ref,control_times_ref); % re-sampled reference profile
            
            % take process measurements at control times of current cycle
            control_times_current_cycle=y.pos4.(['cycle_' num2str(n_rotation-3)]).t_drying((1+...
                p.control_interval/p.drying_sampling_time):(p.control_interval/...
                p.drying_sampling_time):(t/p.drying_sampling_time+1));
            y_dryer_out_Tg=y.pos4.(['cycle_' num2str(n_rotation-3)]).Tg((1+...
                p.control_interval/p.drying_sampling_time):(p.control_interval/...
                p.drying_sampling_time):(t/p.drying_sampling_time+1));
            
            % calculate error
            length_meas=length(y_dryer_out_Tg); % number of measurement at control times of current cycle
            
            if length_meas>2
                errorK=sp.Tg_pos4(length_meas)-y_dryer_out_Tg(end);
                errorK_minus1=sp.Tg_pos4(length_meas-1)-y_dryer_out_Tg(end-1);
                errorK_minus2=sp.Tg_pos4(length_meas-2)-y_dryer_out_Tg(end-2);
                
            elseif length_meas == 2 % second sampling of current cycle: use error at end of previous cycle, too
                errorK=sp.Tg_pos4(length_meas)-y_dryer_out_Tg(end);
                errorK_minus1=sp.Tg_pos4(length_meas-1)-y_dryer_out_Tg(end-1); 
                errorK_minus2=sp.Tg_pos4(end)-y.cont_sign.pos4.Tg(end);  % last error of previous cycle
                
            else % first sampling of current cycle: use errors at end of previous cycle, too
                errorK=sp.Tg_pos4(length_meas)-y_dryer_out_Tg(end);
                errorK_minus1=sp.Tg_pos4(end)-y.cont_sign.pos4.Tg(end); % last error of previous cycle                
                if length(sp.Tg_pos4) == 1 % check if re-sampled reference profile contains more than one measurement
                    % yes: errorK_minus2 taken from two cycles ago
                    errorK_minus2=sp.Tg_pos4(end)-y.cont_sign.pos4.Tg(end-p.control_interval/p.drying_sampling_time-1);
                else
                    %no: errorK_minus2 taken from two cycles ago
                    errorK_minus2=sp.Tg_pos4(end-1)-y.cont_sign.pos4.Tg(end-p.control_interval/p.drying_sampling_time);
                end
                
            end
            
            if t>u.t_rot*.9 % store error in the last 10% of the batch for cycle duration controller
                u.error_t_rot(n_rotation-3) = u.error_t_rot(n_rotation-3) +errorK;
            end
            
        %% MVs update
            if abs(errorK)>.1
                
                % first range: act on pressure drop - MV=dP, CV=Tg,out_drying
                if u.Tinlet_drying <= u_ss.Tinlet_drying   
                    Kc=.0001;
                    tauI=.0001;
                    tauD=0;
                    u.dP_drying=max(3e4,min(u.dP_drying+Kc*(errorK-errorK_minus1+p.control_interval/tauI*errorK+tauD/p.control_interval*...
                            (errorK-2*errorK_minus1+errorK_minus2)), 6e4));
                    u.dP=u.dP_drying;

                end
                
                % second range: act on drying inlet temperature - MV=Tg,in_drying, CV=Tg,out_drying
                if u.dP_drying>=6e4 
                    Kc=.00005;
                    tauI=.05;
                    tauD=0;

                    u.Tinlet_drying=max(u_ss.Tinlet_drying,min(u.Tinlet_drying+Kc*(errorK-errorK_minus1+2/tauI*errorK+tauD/2*...
                    (errorK-2*errorK_minus1+errorK_minus2)),70+273.15));   
                    if u.t_rot > 180
                        u.Tinlet_drying=70+273.15;
                    end
                end
            end
        end
        
    end
end