function simulation_output = run_simulation(u,...
    cryst_output,disturbance_flag,control_flag,total_duration,...
    control_interval,sampling_time)     

%% Check sampling and control times and cycle duration
    if abs(round(1/sampling_time)-1/sampling_time)>0
        error('Sampling time must be a submultiple of 1')
    end
    if abs(round(control_interval)-control_interval)>0
        error('Control interval must be a multiple of 1')
    end
    if abs(round(u.t_rot)-u.t_rot)>0
        error('Cycle duration must be an integer')
    end

%% Initialization           
    % initialize crystallization feed
    load CSD
    cryst_output.x=x;
    cryst_output.CSD=CSD/100/(pi/6)./x.^4; % volume based distribution
    cryst_output.CSD_perc=CSD_perc/100; % volume percentage distribution 
    clear x
    cryst_output.T=295.25;          % (K) - equal to room temperature
    cryst_output_nominal=cryst_output;
    cryst_output.conc_slurry_vector=[];

    % Initialize parameters object
    p=carousel_parameters_class; % create object storing physical properties and operating conditions
    p.integration_interval=1;
    p.control_interval=control_interval;
    p.filtration_sampling_time=sampling_time;
    p.drying_sampling_time=sampling_time;
    p.T_room=cryst_output.T;
    
    % Initialize disturbances vector
    d.V_slurry_dist=1; 
    d.c_slurry_dist=1; 
    d.E_dist=1; 
    d.alpha_dist=1;     

    % initialize disturbances profiles
    d.c_slurry=1+randn(1,1000)*0.02; 
    d.V_slurry=1+randn(1,1000)*0.02;  
    d.E=1+randn(1,1000)*0.02; 
    d.alpha=1+(1-d.E)./d.E.^3;  
    d.hM=1+randn(1,1000)*0.02; 
    d.hT=1+randn(1,1000)*0.02; 
    d.fouling=load('resistances');
    d.resistances=d.fouling.resistances;
    d.stations_working=d.fouling.stations_working;
 
    % Create parameters object    
    p.min_length_discr=3e-4;   % grid spacing for deliquoring and drying
    p = carousel_parameters(p);     
    
    % Initialize object containing manipulated variables profile
    % updated online
    manipulated_vars.t_vector=[]; % process time vector
    manipulated_vars.dP_vector=[];
    manipulated_vars.Tin_drying_vector=[];
    % updated at the end of every cycle
    manipulated_vars.n_cycle_vector=[];
    manipulated_vars.t_rot_vector=[];  % rotation time vector
    manipulated_vars.V_slurry_vector=[]; % fed slurry vector
    
    % time variables initialization
    process_time = 0;
    n_cycle = 0;
    cycle_time = 0;    
   
    %% Sensors object creation
    % measurements
    measurements.t_meas=0;
    measurements.m_filt_WI101=0;
    measurements.P_PI102=1e5;
    measurements.c_slurry_AI101=250;
    measurements.L_cake_LI101=0;
    measurements.V_slurry_LI101=0;            
    measurements.Tg_in_TI101=round(u.Tinlet_drying,1);
    measurements.Tg_out_TI102=round(p.T_room,1);
    measurements.Vdryer_FI101=0;
    
    % noise-free measurements
    measurements_nf.t_meas=0;
    measurements_nf.m_filt_WI101=0;
    measurements_nf.P_PI102=1e5;
    measurements_nf.c_slurry_AI101=250;
    measurements_nf.L_cake_LI101=0;
    measurements_nf.V_slurry_LI101=0;            
    measurements_nf.Tg_in_TI101=round(u.Tinlet_drying,1);
    measurements_nf.Tg_out_TI102=round(p.T_room,1);
    measurements_nf.Vdryer_FI101=0; 
    
    %% other inizializations
    u_nominal=u;
    x=[];   
    y.cake_counter=zeros(1,4);
    y.final_content=[];
    x_estim=[];
    
    %% Simulation
    while process_time <= total_duration 
        
        if cycle_time < u.t_rot && n_cycle > 0   % simulate a step of process operation
            
           % duration of the next simulation step
           simulation_step = min(1, u.t_rot - cycle_time); 
           
           % simulation
           [x,y,measurements,measurements_nf]=carousel_simulator(cycle_time,simulation_step,p,d,u,x,y,measurements,measurements_nf,n_cycle);     
                       
           % update timers
           cycle_time = cycle_time + simulation_step; 
           process_time = process_time + simulation_step;
           
           % call online estimation routines
           x_estim = estimator_online(process_time,cycle_time,...
               p.stations_working,u,measurements,manipulated_vars,x_estim,...
               n_cycle,control_flag,p.filtration_sampling_time,...
               p.control_interval,simulation_step);
           
           % call online control routines and save MVs profiles
           if ceil(cycle_time/p.control_interval)== cycle_time/p.control_interval         
               % other control strategies
               [u,manipulated_vars] = controller_online(process_time,cycle_time,...
               p.stations_working,u,u_nominal,cryst_output_nominal,measurements,manipulated_vars,x_estim,...
               n_cycle,control_flag);
           end
           
        else % rotation                      
           
           % call end of cycle estimation routines
           x_estim = estimator_cycle_switch(process_time,cycle_time,...
               p.stations_working,u,measurements,manipulated_vars,x_estim,...
               n_cycle,control_flag,p.filtration_sampling_time,p.control_interval);
           
           % calculate ethanol content in discharged cake
           y = final_composition(p,x,y,n_cycle);
           
           % update timers
           cycle_time = 0;
           n_cycle = n_cycle+1;
           p.stations_working = d.stations_working(n_cycle,:);
           
           % call end of cycle control routines and save MVs profiles
           [u,u_nominal,manipulated_vars] = controller_cycle_switch(process_time,cycle_time,...
               p.stations_working,u,u_nominal,cryst_output_nominal,measurements,...
               manipulated_vars,x_estim,n_cycle,control_flag);
           
           % call disturbance function
           [cryst_output,d,p]=disturbances(process_time,cryst_output,cryst_output_nominal,p,d,u,n_cycle,disturbance_flag);

           % call switch cycle routine
           [x,y]=switch_cycle(process_time,cryst_output_nominal,p,d,u,x,y,measurements,n_cycle);
       end
    end

%% Prepare output object
if length(manipulated_vars.n_cycle_vector)>4
    d1.resistances=d.resistances(1:manipulated_vars.n_cycle_vector(end),:);
    d1.c_slurry=d.c_slurry(1:manipulated_vars.n_cycle_vector(end));
    d1.V_slurry=d.V_slurry(1:manipulated_vars.n_cycle_vector(end));
    d1.E=d.E(1:manipulated_vars.n_cycle_vector(end));
    d1.alpha=d.alpha(1:manipulated_vars.n_cycle_vector(end));
    d1.hM=d.hM(1:manipulated_vars.n_cycle_vector(end));
    d1.hT=d.hT(1:manipulated_vars.n_cycle_vector(end));

    simulation_output.states=y.states;
    simulation_output.measurements=measurements;
    simulation_output.measurements_nf=measurements_nf;
    simulation_output.disturbances=d1;
    simulation_output.operating_vars=manipulated_vars;
    simulation_output.x_estim=x_estim;
    simulation_output.feed.c_slurry_nominal_vector=cryst_output.conc_slurry_vector;
    simulation_output.cakes_proc_times=y.processing_times;
    simulation_output.final_content=y.final_content;
    simulation_output.active_stations=d.stations_working(1:manipulated_vars.n_cycle_vector(end),:);

    simulation_output.settings.control_mode=control_flag;
    simulation_output.settings.operating_mode=disturbance_flag;
    simulation_output.settings.control_interval=control_interval;
    simulation_output.settings.sampling_time=sampling_time;
    simulation_output.settings.total_duration=total_duration;
    simulation_output.settings.cryst_output_nom=cryst_output_nominal;
    simulation_output.settings.u_nom=u_nominal;
else
    disp('No cakes discharged: increase simulation duration')
    simulation_output='';
end