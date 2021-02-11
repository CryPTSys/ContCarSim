function u = controller(sp,t,u_ss,p,u,y,n_rotation, control_flag)
    
    if control_flag > 0

        %% split range control - MVs=Dp_drying/Tg,out_drying, CV=Tg,out_drying
        if n_rotation>4 % starts from second drying cycle
            % set-point profiles interpolation
            sp.Tg_pos4=interp1(sp.t_ref_Tg,sp.Tg_pos4_ref,(p.control_time*sp.t_rot_ref/u.t_rot):...
                (p.control_time*sp.t_rot_ref/u.t_rot):(u.t_rot*sp.t_rot_ref/u.t_rot));

    %         t_sp=(p.control_time*sp.t_rot_ref/u.t_rot):...
    %             (p.control_time*sp.t_rot_ref/u.t_rot):(u.t_rot*sp.t_rot_ref/u.t_rot);
    %         t_meas=y.pos4.(['cycle_' num2str(n_rotation-3)]).t_drying((1+p.control_time/p.drying_time_step):p.control_time/p.drying_time_step:(t/p.drying_time_step+1));


            y_dryer_out_Tg=y.pos4.(['cycle_' num2str(n_rotation-3)]).Tg((1+p.control_time/p.drying_time_step):p.control_time/p.drying_time_step:(t/p.drying_time_step+1));
            length_meas=length(y_dryer_out_Tg);
            if length_meas>2
                measurement=y_dryer_out_Tg(end-2:end);
                errorK=sp.Tg_pos4(length_meas)-measurement(end);
                errorK_minus1=sp.Tg_pos4(length_meas-1)-measurement(end-1);
                errorK_minus2=sp.Tg_pos4(length_meas-2)-measurement(end-2);
            elseif length_meas == 2
                measurement=y_dryer_out_Tg(end-1:end);
                errorK=sp.Tg_pos4(length_meas)-measurement(end);
                errorK_minus1=sp.Tg_pos4(length_meas-1)-measurement(end-1);
                errorK_minus2=sp.Tg_pos4(end)-y.cont_sign.pos4.Tg(end);
            else
                measurement=y_dryer_out_Tg(end);
                errorK=sp.Tg_pos4(length_meas)-measurement(end);
                errorK_minus1=sp.Tg_pos4(end)-y.cont_sign.pos4.Tg(end);
                if length(sp.Tg_pos4) == 1
                    errorK_minus2=sp.Tg_pos4(end)-y.cont_sign.pos4.Tg(end-p.control_time/p.drying_time_step-1);
                else
                    errorK_minus2=sp.Tg_pos4(end-1)-y.cont_sign.pos4.Tg(end-p.control_time/p.drying_time_step);
                end
            end
            if t>u.t_rot*.9 %&& t>u.t_rot*.5
                u.error_t_rot(n_rotation-3) = u.error_t_rot(n_rotation-3) +errorK;
            end
            if abs(errorK)>.1

                if u.Tinlet_drying <= u_ss.Tinlet_drying   
                    % PID 
                    Kc=.0001;
                    tauI=.0001;
                    tauD=0;
                    u.dP_drying=max(3e4,min(u.dP_drying+Kc*(errorK-errorK_minus1+p.control_time/tauI*errorK+tauD/p.control_time*...
                            (errorK-2*errorK_minus1+errorK_minus2)),6e4));

                end
                if u.dP_drying>=6e4 % MV=Tg,in_drying, CV=Tg,out_drying
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