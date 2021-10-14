function x_estim = estimator_cycle_switch(x_estim,u,measurements,control_flag)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs
    %
    % x_estim           =   object containing states and parameters estimated by estimator_online.m and estimator_cycle_end
    %                       Fields follow the structure defined in run_carousel.m
    % u                 =   vector of manipulated variables during previous control interval
    %                       Fields of u:
    %                       - u.t_rot=cycle duration (s)
    %                       - u.V_slurry=fed slurry volume (m3)
    %                       - u.dP=pressure drop Stations 1-4 (Pa)
    %                       - u.Tinlet_drying=drying gas temperature Station 5 (K)                      
    % measurements      =   object of process measurements since process onset 
    %                       with the sampling time that has been set in run_carousel.m
    %                       Fields of measurements:
    %                       - measurements.t = vector of sampling times (s)
    %                       - measurements.m_filt_WI101 = vector of filtrate mass measured by WI101 (kg)
    %                       - measurements.P_PI102 = vector of pressure measured by PI102 (Pa)
    %                       - measurements.c_slurryAI101 = vector of slurry concentration measured by AI101 (kg/m3)
    %                       - measurements.L_cake_LI101 = vector of height of cakes in Station 1 measured by LI101 (m)
    %                       - measurements.V_slurry_LI101 = vector of slurry volume in Station 1 measured by LI101 (m)
    %                       - measurements.Tg_top_TI101 = vector of temperatures of drying gas measured by TI101 (K) - inlet
    %                       - measurements.Tg_bot_TI102 = vector of temperatures of drying gas measured by TI102 (K) - outlet                      
    %                       - measurements.Vdryer_FI101 = vector of drying gas flowrate measured by FI101 (m3/s)
    % control_flag      =   scalar defined in run_carousel.m
    %
    % Outputs           
    % x_estim           =   Updated after call to this function 
    %                       Object containing states and parameters
    %                       estimated by estimator_online.m and estimator_cycle_switch
    %                       Fields follow the structure defined in run_carousel.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
end

