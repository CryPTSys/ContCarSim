function [u,operating_vars] = controller_online(process_time,cycle_time,...
    stations_working,u,u_nominal,cryst_output_nominal,measurements,operating_vars,x_estim,n_cycle,control_mode)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs
    %
    % process_time      =   timer started at process onset (s)
    % cycle_time        =   timer re-started at every carousel rotation (s)
    % stations_working     =   vector [1x4] - for i=1:4:
    %                       - stations_working(i)=1 if during current cycle station i is processing material;
    %                       - stations_working(i)=0 if during current cycle station i is empty; 
    % u                 =   vector of manipulated variables during previous control interval
    %                       Fields of u:
    %                       - u.t_rot=cycle duration set-point (s)    MUST BE AN INTEGER
    %                       - u.V_slurry=fed slurry volume set-point (m3)
    %                       - u.dP=pressure drop Stations 1-4 set point (Pa)
    %                       - u.Tinlet_drying=drying gas temperature Station 5 set-point (K)                       
    % u_nominal         =   nominal value of manipulated variables, as set in run_carousel.m
    %                       Same fields of u
    % cryst_output_nominal  = object containing nominal (i.e., initial) feed
    %                        conditions. Fields:
    %                        - cryst_output_nominal.conc_slurry= nominal
    %                                 slurry concentration in feed (kg/m3)
    %                        - cryst_output_nominal.x = Crystal size
    %                                 distribution – particles diameters (m)
    %                        - cryst_output_nominal.CSD - Volumetric crystal size distribution – percentage
    %                        - cryst_output_nominal.CSD_perc - Volumetric crystal size distribution – percentage
    %                        - cryst_output_nominal.T - slurry temperature (= room temperature) (K)  
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
    %                       - operating_vars.n_cycle_vector = list of number of initialized carousel cycles
    %                       - operating_vars.t_rot_vector = u.t_rot time profile [1 x length(operating_vars.n_cycle_vector)]
    %                       - operating_vars.V_slurry_vector = u.V_slurry time profile [1 x length(operating_vars.n_cycle_vector)]
    % x_estim           =   object containing states and parameters estimated by estimator_online.m and estimator_cycle_switch
    %                       Fields follow the structure defined in run_carousel.m
    % n_cycle           =   cycle counter - number of current cycle
    % control_mode      =   scalar defined in run_carousel.m
    %
    % Outputs           
    % u                 =   vector of manipulated variables for follwing control interval
    %                       Fields of u:
    %                       - u.t_rot=cycle duration set-point (s)    MUST BE AN INTEGER
    %                       - u.V_slurry=fed slurry volume set-point (m3)
    %                       - u.dP=pressure drop Stations 1-4 set point (Pa)
    %                       - u.Tinlet_drying=drying gas temperature Station 5 set-point (K)   
    %         -------->     Fields not updated during call to this function
    %                       will retain the value set for the previous control interval
    % operating_vars  =   object storing the profiles of the manipulated variables (automatically updated)
    %                       Fields of operating_vars:
    %                       - operating_vars.t_vector = control times vector
    %                       - operating_vars.dP_vector = u.dP time profile [1 x length(operating_vars.t_vector)]
    %                       - operating_vars.Tin_drying_vector = u.Tinlet_drying time profile [1 x length(operating_vars.t_vector)]
    %                       - operating_vars.n_cycle_vector = list of number of initialized carousel cycles
    %                       - operating_vars.t_rot_vector = u.t_rot time profile [1 x length(operating_vars.n_cycle_vector)]
    %                       - operating_vars.V_slurry_vector = u.V_slurry time profile [1 x length(operating_vars.n_cycle_vector)]
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %% Paper simulator
    % note that u.t_rot must always be an integer!
    
    if control_mode == 1 && stations_working(4)==1 % Sample closed-loop control - works when Station 4 is working
        if measurements.Tg_out_TI102(end)<295.3 || ...
            measurements.Tg_out_TI102(end)<=measurements.Tg_out_TI102(end-2)               
            u.t_rot=cycle_time+1;
        else % trigger cycle end
            u.t_rot=cycle_time;
        end
    elseif control_mode == 0 || stations_working(4)==0 % open-loop - or if Station 4 is empty
        u.t_rot=u_nominal.t_rot;
    end   
   %% do not modify part below
   % Store manipulated variables profile
   operating_vars.t_vector=[operating_vars.t_vector process_time];     
   operating_vars.dP_vector=[operating_vars.dP_vector u.dP];
   operating_vars.Tin_drying_vector=[operating_vars.Tin_drying_vector u.Tinlet_drying];
   
end