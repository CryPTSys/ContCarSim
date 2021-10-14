function [u,manipulated_vars] = controller_cycle_switch(process_time,batch_time,...
               ports_working,u,u_nominal,measurements,manipulated_vars,x_estim,...
               n_cycle,control_flag)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs
    %
    % process_time      =   timer started at process onset (s)
    % batch_time        =   timer re-started at every carousel rotation (s)
    % ports_working     =   vector [1x4] - for i=1:4:
    %                       - ports_working(i)=1 if during current cycle port i is processing material;
    %                       - ports_working(i)=0 if during current cycle port i is empty; 
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
    %                       - measurements.Tg_top_TI101 = vector of temperatures of drying gas measured by TI101 (K) - inlet
    %                       - measurements.Tg_bot_TI102 = vector of temperatures of drying gas measured by TI102 (K) - outlet                      
    %                       - measurements.Vdryer_FI101 = vector of drying gas flowrate measured by FI101 (m3/s)
    % manipulated_vars  =   object storing the profiles of the manipulated variables (automatically updated)
    %                       Fields of manipulated_vars:
    %                       - manipulated_vars.t_vector = control times vector
    %                       - manipulated_vars.dP_vector = u.dP time profile [1 x length(manipulated_vars.t_vector)]
    %                       - manipulated_vars.Tin_drying_vector = u.Tinlet_drying time profile [1 x length(manipulated_vars.t_vector)]
    %                       - manipulated_vars.n_cycle_vector = list of number of completed carousel cycles
    %                       - manipulated_vars.t_rot_vector = u.t_rot time profile [1 x length(manipulated_vars.n_cycle_vector)]
    %                       - manipulated_vars.V_slurry_vector = u.V_slurry time profile [1 x length(manipulated_vars.n_cycle_vector)]
    % x_estim           =   object containing states and parameters estimated by estimator_online.m and estimator_cycle_switch
    %                       Fields follow the structure defined in run_carousel.m
    % n_cycle           =   cycle counter
    % control_flag      =   scalar defined in run_carousel.m
    %
    % Outputs           
    % u                 =   vector of manipulated variables for follwing control interval
    %                       Fields of u:
    %                       - u.t_rot=cycle duration (s)
    %                       - u.V_slurry=fed slurry volume (m3)
    %                       - u.dP=pressure drop Stations 1-4 (Pa)
    %                       - u.Tinlet_drying=drying gas temperature Station 5 (K)   
    %         -------->     Fields not updated during call to this function
    %                       will retain the value set for the previous control interval
    % manipulated_vars  =   object storing the profiles of the manipulated variables (automatically updated)
    %                       Fields of manipulated_vars:
    %                       - manipulated_vars.t_vector = control times vector
    %                       - manipulated_vars.dP_vector = u.dP time profile [1 x length(manipulated_vars.t_vector)]
    %                       - manipulated_vars.Tin_drying_vector = u.Tinlet_drying time profile [1 x length(manipulated_vars.t_vector)]
    %                       - manipulated_vars.n_cycle_vector = list of number of completed carousel cycles
    %                       - manipulated_vars.t_rot_vector = u.t_rot time profile [1 x length(manipulated_vars.n_cycle_vector)]
    %                       - manipulated_vars.V_slurry_vector = u.V_slurry time profile [1 x length(manipulated_vars.n_cycle_vector)]
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    
    %% Layers 0-1
    if ports_working(1)==0
       u.V_slurry=0;
    else
       u.V_slurry=u_nominal.V_slurry;
    end

    %% Store manipulated variables profile

    manipulated_vars.n_cycle_vector=[manipulated_vars.n_cycle_vector n_cycle];         
    manipulated_vars.V_slurry_vector=[manipulated_vars.V_slurry_vector u.V_slurry];
    if n_cycle > 1
        manipulated_vars.t_rot_vector=[manipulated_vars.t_rot_vector u.t_rot];
    end
end