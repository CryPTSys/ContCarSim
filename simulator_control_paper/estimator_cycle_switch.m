function x_estim = estimator_cycle_switch(process_time,cycle_time,...
               stations_working,u,measurements,operating_vars,x_estim,...
               n_cycle,control_mode, sampling_time, control_interval)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function estimating parameters and states at every carousel rotation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs
    %
    % process_time      =   timer started at process onset (s)
    % cycle_time        =   timer re-started at every carousel rotation (s)
    % stations_working     =   vector [1x4] - for i=1:4:
    %                       - stations_working(i)=1 if station i was processing material during cycle that just finished ;
    %                       - stations_working(i)=0 if station i was empty during cycle that just finished; 
    % u                 =   vector of manipulated variables during previous control interval
    %                       Fields of u:
    %                       - u.t_rot=cycle duration (s)
    %                       - u.V_slurry=fed slurry volume (m3)
    %                       - u.dP=pressure drop Stations 1-4 (Pa)
    %                       - u.Tinlet_drying=drying gas temperature Station 5 (K)                      
    % u_nominal         =   nominal value of manipulated variables, as set in run_carousel.m
    %                       Same fields of u
    % measurements      =   object of process measurements since process onset 
    %                       with the sampling time that has been set in run_carousel.m
    %                       Fields of measurements:
    %                       - measurements.t_meas = vector of sampling times (s)
    %                       - measurements.m_filt_WI101 = vector of filtrate mass measured by WI101 (kg)
    %                       - measurements.P_PI102 = vector of pressure measured by PI102 (Pa)
    %                       - measurements.c_slurryAI101 = vector of slurry concentration measured by AI101 (kg/m3)
    %                       - measurements.L_cake_LI101 = vector of height of cakes in Station 1 measured by LI101 (m)
    %                       - measurements.V_slurry_LI101 = vector of slurry volume in Station 1 measured by LI101 (m)
    %                       - measurements.Tg_in_TI101 = vector of temperatures of drying gas measured by TI101 (K) - inlet
    %                       - measurements.Tg_out_TI102 = vector of temperatures of drying gas measured by TI102 (K) - outlet                      
    %                       - measurements.Vdryer_FI101 = vector of drying gas flowrate measured by FI101 (m3/s)
    % operating_vars  =   object storing the profiles of the manipulated variables (automatically updated)
    %                       Fields of operating_vars:
    %                       - operating_vars.t_vector = control times vector
    %                       - operating_vars.dP_vector = u.dP time profile [1 x length(operating_vars.t_vector)]
    %                       - operating_vars.Tin_drying_vector = u.Tinlet_drying time profile [1 x length(operating_vars.t_vector)]
    %                       - operating_vars.n_cycle_vector = list of number of completed carousel cycles
    %                       - operating_vars.t_rot_vector = u.t_rot time profile [1 x length(operating_vars.n_cycle_vector)]
    %                       - operating_vars.V_slurry_vector = u.V_slurry time profile [1 x length(operating_vars.n_cycle_vector)]
    % x_estim           =   object containing states and parameters estimated by estimator_online.m and estimator_cycle_switch
    %                       Fields follow the structure defined in run_carousel.m
    % n_cycle           =   cycle counter - number of cycle that has just finished
    % control_mode      =   scalar defined in run_carousel.m
    % sampling_time     =   sampling time (s)
    % control_interval  =   control_interval (s)
    %
    % Outputs           
    % x_estim           =   Updated after call to this function 
    %                       Object containing states and parameters
    %                       estimated by estimator_online.m and estimator_cycle_switch
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
end
