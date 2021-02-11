function u = controller_slow(sp,t,u_ss,p,u,y,n_rotation, control_flag)
    if control_flag > 0
        if n_rotation>5 
    %         % set-point profiles interpolation
    %         sp.Tg_pos5=interp1(sp.t_ref_Tg,sp.Tg_pos5_ref,(p.control_time*sp.t_rot_ref/u.t_rot):...
    %             (p.control_time*sp.t_rot_ref/u.t_rot):(u.t_rot*sp.t_rot_ref/u.t_rot));
    % 
    %         y_dryer_out_Tg=y.pos5.(['cycle_' num2str(n_rotation-4)]).Tg((1+p.control_time/p.drying_time_step):p.control_time/p.drying_time_step:(t/p.drying_time_step+1));
    %         length_meas=length(y_dryer_out_Tg);
    %         if length_meas>2
    %             measurement=y_dryer_out_Tg(end-2:end);
    %             errorK=sp.Tg_pos5(length_meas)-measurement(end);
    %             errorK_minus1=sp.Tg_pos5(length_meas-1)-measurement(end-1);
    %             errorK_minus2=sp.Tg_pos5(length_meas-2)-measurement(end-2);
    %         elseif length_meas == 2
    %             measurement=y_dryer_out_Tg(end-1:end);
    %             errorK=sp.Tg_pos5(length_meas)-measurement(end);
    %             errorK_minus1=sp.Tg_pos5(length_meas-1)-measurement(end-1);
    %             errorK_minus2=sp.Tg_pos5(end)-y.cont_sign.pos5.Tg(end);
    %         else
    %             measurement=y_dryer_out_Tg(end);
    %             errorK=sp.Tg_pos5(length_meas)-measurement(end);
    %             errorK_minus1=sp.Tg_pos5(end)-y.cont_sign.pos5.Tg(end);
    %             if length(sp.Tg_pos5) == 1
    %                 errorK_minus2=sp.Tg_pos5(end)-y.cont_sign.pos5.Tg(end-p.control_time/p.drying_time_step-1);
    %             else
    %                 errorK_minus2=sp.Tg_pos5(end-1)-y.cont_sign.pos5.Tg(end-p.control_time/p.drying_time_step);
    %             end
    %         end

    %             if abs(u.error_t_rot(n_rotation-3))>.1
                    errorK=u.error_t_rot(n_rotation-3);
                if n_rotation>5
                    errorK_minus1=u.error_t_rot(n_rotation-4);
                else
                    errorK_minus1=0;
                end
                if n_rotation>6
                    errorK_minus2=u.error_t_rot(n_rotation-5);
                else
                    errorK_minus2=0;
                end

                if u.Tinlet_drying >= 70+273.15 
                    % springer chapter
    %                 Kc=5e-4;
    %                 tauI=.008; % 0.005
    %                 if errorK<0
    %                     tauI=.0015; % even lower if you want to go back to the original rotation time
    %                 end
                    % escape v2
                    Kc=5e-4;
                    tauI=.008; % 0.005
                    if errorK<0
                        tauI=.001; % even lower if you want to go back to the original rotation time
                    end

                    tauD=0;
                    u.t_rot=round(max(1,u.t_rot+Kc*(errorK-errorK_minus1+2/tauI*errorK)));

                elseif u.dP_drying<=3e4
                    Kc=2.8e-4;
                    tauI=.001;
                    u.t_rot=round(max(1,u.t_rot+Kc*(errorK-errorK_minus1+2/tauI*errorK)));
                end

    %         p.control_time=u.t_rot;    
            %% ratio controller
    %         u.flowrate_slurry=u_ss.V_slurry/u.t_rot;%.025e-6;%1.7e-7; % m3/s

            u.error_t_rot(n_rotation+1)=0;
        end
    end
end