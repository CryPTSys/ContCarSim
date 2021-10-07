function [u,controller_output,x_estim] = controller_cycle_end(t,u,y,controller_output,x_estim,n_cycle,control_flag)
   


   %% store manipulated variables profile
   controller_output.n_cycle_vector=[controller_output.n_cycle_vector n_cycle+1];     
   controller_output.t_rot_vector=[controller_output.t_rot_vector u.t_rot];
   controller_output.V_slurry_vector=[controller_output.V_slurry_vector u.V_slurry];
     
end