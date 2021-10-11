function [u,manipulated_vars,x_estim] = controller_online(t,ports_working,u,u_nominal,measurements,manipulated_vars,x_estim,n_cycle,control_flag)
        
    if control_flag==0 % Layer 0
        u.t_rot=u_nominal.t_rot;   
    elseif control_flag == 1 && sum(ports_working==4)==1 % Layer 1
        if measurements.Tg_bot_TI102(end)<18.7+273.15 || ...
            measurements.Tg_bot_TI102(end)<=measurements.Tg_bot_TI102(end-2)               
            u.t_rot=t+1;%p.control_interval;
        else % trigger cycle end
            u.t_rot=t;
        end
    else
        u.t_rot=u_nominal.t_rot;
    end
    
    
   %% store manipulated variables profile
   manipulated_vars.t_vector=[manipulated_vars.t_vector t];     
   manipulated_vars.dP_vector=[manipulated_vars.dP_vector u.dP_drying];
   manipulated_vars.Tin_drying_vector=[manipulated_vars.Tin_drying_vector u.Tinlet_drying];
end