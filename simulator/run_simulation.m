function simulation_output = run_simulation(u,...
    cryst_output,disturbance_scenario,control_mode,total_duration,...
    control_interval,sampling_interval,inter_cycle_Dt,mesh_clean_Dt)     

%% Check sampling and control times and cycle duration
    if abs(round(1/sampling_interval)-1/sampling_interval)>0
        error('Sampling time must be a submultiple of 1')
    end
    if abs(round(control_interval)-control_interval)>0
        error('Control interval must be a multiple of 1')
    end
    if abs(round(u.t_cycle)-u.t_cycle)>0
        error('Cycle duration must be an integer')
    end
    if u.V_slurry>10e-6
        error('Loaded slurry volume exceeds maximum capacity')
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
    p.filtration_sampling_interval=sampling_interval;
    p.drying_sampling_interval=sampling_interval;
    p.T_room=cryst_output.T;
    
    % Initialize parametric disturbances 
    % will be reassigned at the onset of every cycle
    d.V_slurry_dist=1; 
    d.c_slurry_dist=1; 
    d.E_dist=1; 
    d.alpha_dist=1;     

    % Create disturbances profiles
    % values that disturbances will assume cycle after cycle
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
    
    % Initialize object containing operating variables profile
    
    % operating variables stored at every control interval
    operating_vars.t_vector=[]; % process time vector
    operating_vars.P_compr_vector=[];
    operating_vars.Tin_drying_vector=[];
    
    % operating variables stored only at the end of a cycle
    operating_vars.n_cycle_vector=[];
    operating_vars.t_cycle_vector=[];  % rotation time vector
    operating_vars.V_slurry_vector=[]; % fed slurry vector
    
    % time variables initialization
    process_time = 0;
    n_cycle = 0;
    cycle_time = 0;    
   
    %% Sensors object creation
    % initialized to the value they assume at process onset
    % measurements 
    measurements.t_meas=0;
    measurements.m_filt_WI101=0;
    measurements.P_PI101=u.P_compr;
    measurements.P_PI102=0;
    measurements.c_slurry_AI101=cryst_output.conc_slurry;
    measurements.L_cake_LI101=0;
    measurements.V_slurry_LI101=0;            
    measurements.Tg_in_TI101=round(u.Tinlet_drying,1);
    measurements.Tg_out_TI102=round(p.T_room,1);
    measurements.Vdryer_FI101=0;
    
    % noise-free measurements
    measurements_nf.t_meas=0;
    measurements_nf.m_filt_WI101=0;
    measurements_nf.P_PI101=u.P_compr;
    measurements_nf.P_PI102=0;
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
        
        if cycle_time < u.t_cycle && n_cycle > 0   % simulate a step of process operation
           if process_time < total_duration
               % duration of the next simulation step
               simulation_step = min(1, u.t_cycle - cycle_time); 
               
               % simulation
               [x,y,measurements,measurements_nf]=carousel_simulator(cycle_time,simulation_step,p,d,u,x,y,measurements,measurements_nf,n_cycle);     
                           
               % update timers
               cycle_time = cycle_time + simulation_step; 
               process_time = process_time + simulation_step;
               
               % call online estimation routines
               x_estim = estimator_online(process_time,cycle_time,...
                   p.stations_working,u,measurements,operating_vars,x_estim,...
                   n_cycle,control_mode,p.filtration_sampling_interval,...
                   p.control_interval,simulation_step);
    
               % call online control routines and save MVs profiles
               if ceil(cycle_time/p.control_interval)== cycle_time/p.control_interval         
                   [u,operating_vars] = controller_online(process_time,cycle_time,...
                   p.stations_working,u,u_nominal,cryst_output_nominal,measurements,operating_vars,x_estim,...
                   n_cycle,control_mode);
               end
           else
               process_time = process_time+1;
           end

           
        else % rotation                      
           
           % call end of cycle estimation routines
           x_estim = estimator_cycle_switch(process_time,cycle_time,...
               p.stations_working,u,measurements,operating_vars,x_estim,...
               n_cycle,control_mode,p.filtration_sampling_interval,p.control_interval);

           % calculate ethanol content in discharged cake
           y = final_composition(p,x,y,n_cycle);
           
           % update timers
           cycle_time = 0;
           n_cycle = n_cycle+1;
           p.stations_working = d.stations_working(n_cycle,:);

           % call end of cycle control routines and save MVs profiles
           [u,u_nominal,operating_vars] = controller_cycle_switch(process_time,cycle_time,...
               p.stations_working,u,u_nominal,cryst_output_nominal,measurements,...
               operating_vars,x_estim,n_cycle,control_mode);

           % call disturbance function
           [cryst_output,d,p]=disturbances(process_time,cryst_output,cryst_output_nominal,p,d,u,n_cycle,disturbance_scenario);

           if process_time > 0
               if sum(p.stations_working == [1 0 0 0])==4
                   process_time=process_time+mesh_clean_Dt;
               end
               process_time=process_time+inter_cycle_Dt;               
           end

           % call switch cycle routine
           [x,y,measurements,measurements_nf]=switch_cycle(process_time,...
               cryst_output_nominal,p,d,u,x,y,measurements,measurements_nf,...
               operating_vars,n_cycle);

       end
    end
    
    %% Prepare output object
    if cycle_time<=1
        operating_vars.n_cycle_vector(end)=[];
        operating_vars.V_slurry_vector(end)=[];
        cryst_output.conc_slurry_vector(end)=[];
    else
        if cycle_time < u.t_cycle 
            disp('Note: Final cycle not finished: outputs related to final cycle not included in simulation_output.m') %operating_vars.t_cycle_vector(end+1)=cycle_time;
        end
    end

    if length(operating_vars.n_cycle_vector)>4
        d1.resistances=d.resistances(1:operating_vars.n_cycle_vector(end),:);
        d1.c_slurry=d.c_slurry(1:operating_vars.n_cycle_vector(end));
        d1.V_slurry=d.V_slurry(1:operating_vars.n_cycle_vector(end));
        d1.E=d.E(1:operating_vars.n_cycle_vector(end));
        d1.alpha=d.alpha(1:operating_vars.n_cycle_vector(end));
        d1.hM=d.hM(1:operating_vars.n_cycle_vector(end));
        d1.hT=d.hT(1:operating_vars.n_cycle_vector(end));

        simulation_output.states=y.states;
        simulation_output.measurements=measurements;
        simulation_output.measurements_nf=measurements_nf;
        simulation_output.disturbances=d1;
        simulation_output.operating_vars=operating_vars;
        simulation_output.x_estim=x_estim;
        simulation_output.feed.c_slurry_nom_vector=cryst_output.conc_slurry_vector;
        simulation_output.cakes_proc_times=y.processing_times;
        simulation_output.final_content=y.final_content;
        simulation_output.active_stations=d.stations_working(1:operating_vars.n_cycle_vector(end),:);
        
        % throughput calculation
        slurry_volumes=simulation_output.operating_vars.V_slurry_vector.*simulation_output.disturbances.V_slurry;
        null_volumes=(slurry_volumes==0);
        slurry_volumes(null_volumes)=[];
        slurry_volumes=slurry_volumes(1:length(y.final_content));
        slurry_concs=simulation_output.feed.c_slurry_nom_vector.*simulation_output.disturbances.c_slurry;
        slurry_concs(null_volumes)=[];
        slurry_concs=slurry_concs(1:length(y.final_content));
        cakes_mass=slurry_concs.*slurry_volumes;              
        acceptable_cakes=y.final_content<=0.005;
        total_production=sum(cakes_mass(acceptable_cakes)); % kg
        simulation_output.total_production=total_production;
        
        % include input setting in output object
        simulation_output.settings.control_mode=control_mode;
        simulation_output.settings.disturbance_scenario=disturbance_scenario;
        simulation_output.settings.control_interval=control_interval;
        simulation_output.settings.sampling_interval=sampling_interval;
        simulation_output.settings.total_duration=total_duration;
        simulation_output.settings.cryst_output_nom.conc_slurry=cryst_output_nominal.conc_slurry;
        simulation_output.settings.cryst_output_nom.x=cryst_output_nominal.x;
        simulation_output.settings.cryst_output_nom.CSD_perc=cryst_output_nominal.CSD_perc;
        simulation_output.settings.cryst_output_nom.T=cryst_output_nominal.T;
        simulation_output.settings.u_nom=u_nominal;
        simulation_output.settings.inter_cycle_Dt=inter_cycle_Dt;
        simulation_output.settings.mesh_clean_Dt=mesh_clean_Dt;

    else
        warning('No cakes discharged: increase simulation duration')
        simulation_output='';
    end