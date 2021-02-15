function u = controller_cycle_duration(~,~,~,~,u,~,n_cycle, control_flag)
    if control_flag > 0
        
        if n_cycle>4  % start updating cycle duration only after first cake discharge
                errorK=u.error_t_rot(n_cycle-3); % use error of last dried batch
                
                % retrieve error of previously dried batches for PID law
                if n_cycle>5
                    errorK_minus1=u.error_t_rot(n_cycle-4);
                else
                    errorK_minus1=0;
                end
                if n_cycle>6
                    errorK_minus2=u.error_t_rot(n_cycle-5);
                else
                    errorK_minus2=0;
                end
                
                % if dP and drying inlet temperature are saturated, update cycle duration
                if u.Tinlet_drying >= 70+273.15 % right part of splitted ranges
                    Kc=5e-4;
                    tauI=.008;
                    u.t_rot=round(max(1,u.t_rot+Kc*(errorK-errorK_minus1+2/tauI*errorK)));

                elseif u.dP_drying<=3e4 % left part of splitted ranges
                    Kc=2.8e-4;
                    tauI=.001;
                    u.t_rot=round(max(1,u.t_rot+Kc*(errorK-errorK_minus1+2/tauI*errorK)));
                end
                
            % add to error vector the element for the next drying batch
            u.error_t_rot(n_cycle+1)=0;
        end
    end
end