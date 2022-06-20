function [u,u_nominal,operating_vars] = controller_cycle_switch_Layer2(process_time,cycle_time,...
               stations_working,u,u_nominal,cryst_output_nominal,measurements,operating_vars,x_estim,...
               n_cycle,control_mode)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs
    %
    % process_time      =   timer started at process onset (s)
    % cycle_time        =   timer re-started at every carousel rotation (s)
    % stations_working     =   vector [1x4] - for i=1:4:
    %                       - stations_working(i)=1 if station i was processing material during cycle that just finished ;
    %                       - stations_working(i)=0 if station i was empty during cycle that just finished; 
    % u                 =   vector of manipulated variables during previous control interval
    %                       Fields of u:
    %                       - u.t_cycle=cycle duration set-point (s)    MUST BE AN INTEGER
    %                       - u.V_slurry=fed slurry volume set-point (m3)
    %                       - u.P_compr= gauge pressure provided by compressor P101 (Pa)
    %                       - u.Tinlet_drying=drying gas temperature Station 5 set-point (K)                   
    % u_nominal         =   nominal value of manipulated variables, as set in run_carousel.m
    %                       or updated later. Same fields of u
    % cryst_output_nominal = object containing nominal (i.e., initial) feed
    %                        conditions. Fields
    %                       - cryst_output_nominal
    % measurements      =   object of process measurements since process onset 
    %                       with the sampling interval that has been set in run_carousel.m
    %                       Fields of measurements:
    %                       - measurements.t_meas = vector of sampling intervals (s)
    %                       - measurements.m_filt_WI101 = vector of filtrate mass measured by WI101 (kg)
    %                       - measurements.P_PI101 = vector of pressure measured by PI101 - gauge (Pa)
    %                       - measurements.P_PI102 = vector of pressure measured by PI102 - gauge (Pa)
    %                       - measurements.c_slurry_AI101 = vector of slurry concentration measured by AI101 (kg/m3)
    %                       - measurements.L_cake_LI101 = vector of height of cakes in Station 1 measured by LI101 (m)
    %                       - measurements.V_slurry_LI101 = vector of slurry volume in Station 1 measured by LI101 (m)
    %                       - measurements.Tg_in_TI101 = vector of temperatures of drying gas measured by TI101 (K) - inlet
    %                       - measurements.Tg_out_TI102 = vector of temperatures of drying gas measured by TI102 (K) - outlet                      
    %                       - measurements.Vdryer_FI101 = vector of drying gas flowrate measured by FI101 (m3/s)
    % operating_vars  =   object storing the profiles of the manipulated variables (automatically updated)
    %                       Fields of operating_vars:
    %                       - operating_vars.t_vector = control times vector
    %                       - operating_vars.P_compr_vector = u.P_compr time profile [1 x length(operating_vars.t_vector)]
    %                       - operating_vars.Tin_drying_vector = u.Tinlet_drying time profile [1 x length(operating_vars.t_vector)]
    %                       - operating_vars.n_cycle_vector = list of number of initialized carousel cycles
    %                       - operating_vars.t_cycle_vector = u.t_cycle time profile [1 x length(operating_vars.n_cycle_vector)]
    %                       - operating_vars.V_slurry_vector = u.V_slurry time profile [1 x length(operating_vars.n_cycle_vector)]
    % x_estim           =   object containing states and parameters estimated by estimator_online.m and estimator_cycle_switch
    %                       Fields follow the structure defined in run_carousel.m
    % n_cycle           =   cycle counter - number of cycle that has just finished
    % control_mode      =   scalar defined in run_carousel.m
    %
    % Outputs           
    % u                 =   vector of manipulated variables for following control interval
    %                       Fields of u:
    %                       - u.t_cycle=cycle duration set-point (s)    MUST BE AN INTEGER
    %                       - u.V_slurry=fed slurry volume set-point (m3)
    %                       - u.P_compr= gauge pressure provided by compressor P101 (Pa)
    %                       - u.Tinlet_drying=drying gas temperature Station 4 set-point (K)   
    %         -------->     Fields not updated during call to this function
    %                       will retain the value set for the previous control interval
    % u_nominal         =   nominal value of manipulated variables, as set in run_carousel.m
    %                       or updated later. Same fields of u
    % operating_vars  =   object storing the profiles of the manipulated variables (automatically updated)
    %                       Fields of operating_vars:
    %                       - operating_vars.t_vector = control times vector
    %                       - operating_vars.P_compr_vector = u.P_compr time profile [1 x length(operating_vars.t_vector)]
    %                       - operating_vars.Tin_drying_vector = u.Tinlet_drying time profile [1 x length(operating_vars.t_vector)]
    %                       - operating_vars.n_cycle_vector = list of number of initialized carousel cycles
    %                       - operating_vars.t_cycle_vector = u.t_cycle time profile [1 x length(operating_vars.n_cycle_vector)]
    %                       - operating_vars.V_slurry_vector = u.V_slurry time profile [1 x length(operating_vars.n_cycle_vector)]
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
    

    %% Layers 0-1-2
    if stations_working(1)==0
       u.V_slurry=0;
    elseif control_mode==3.3 
       opt_conditions=run_rto_robust(stations_working,x_estim,u,u_nominal,...
           cryst_output_nominal,measurements,n_cycle);
       u.V_slurry=opt_conditions(2);               
       % for the first three cycles after mesh cleaning, the rotation time
       % is kept to the optimal value found right after mesh cleaning
       if stations_working(2)==0 
           u_nominal.t_cycle=round(opt_conditions(1));
       end
    else
       u.V_slurry=u_nominal.V_slurry;
    end
    
     %% do not modify part below
    if stations_working(1)==0 % if Station 1 is empty at next cycle
       u.V_slurry=0;          % then no slurry loaded in Station 1
    end

    % Store manipulated variables profile
    operating_vars.n_cycle_vector=[operating_vars.n_cycle_vector n_cycle];         
    operating_vars.V_slurry_vector=[operating_vars.V_slurry_vector u.V_slurry];
    if n_cycle > 1
        operating_vars.t_cycle_vector=[operating_vars.t_cycle_vector u.t_cycle];
    end
    
    if round(u.t_cycle)-u.t_cycle>1e-6
       error('u.t_cycle must be an integer!')
    end
end