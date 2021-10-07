function [u,controller_output,x_estim] = controller_online(t,u,y,controller_output,x_estim,n_cycle,control_flag)
    
    
   %% store manipulated variables profile
   controller_output.t_vector=[controller_output.t_vector t];     
   controller_output.dP_vector=[controller_output.dP_vector u.dP_drying];
   controller_output.Tin_drying_vector=[controller_output.Tin_drying_vector u.Tinlet_drying];
end