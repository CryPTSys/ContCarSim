function [u,manipulated_vars,x_estim] = controller_cycle_switch(t,ports_working,u,u_nominal,measurements,manipulated_vars,x_estim,n_cycle,control_flag)
   
   if sum(ports_working==1)==0
       u.V_slurry=0;
   else
       u.V_slurry=u_nominal.V_slurry;
   end
   
   % store manipulated variables profile
   if n_cycle > 1
      manipulated_vars.n_batch_vector=[manipulated_vars.n_batch_vector n_cycle-1];     
      manipulated_vars.t_rot_vector=[manipulated_vars.t_rot_vector u.t_rot];
      manipulated_vars.V_slurry_vector=[manipulated_vars.V_slurry_vector u.V_slurry];
   end
end