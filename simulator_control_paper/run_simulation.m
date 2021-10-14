function simulation_output = run_simulation(u,...
    cryst_output,disturbance_flag,control_flag,total_duration,...
    manipulated_vars,x_estim,control_interval,sampling_time)     

%% Initialization           
    % Initialize parameters object
    p=carousel_parameters_class_control; % create object storing physical properties and operating conditions
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
    d.hM=1+randn(1,1000)*0.02; 
    d.hT=1+randn(1,1000)*0.02; 
    d.fouling=load('resistances');
    d.resistances=d.fouling.resistances;
    d.ports_working=d.fouling.ports_working;
 
    % initialize crystallization feed
    load CSD
    cryst_output.x=x;
    cryst_output.CSD=CSD/100/(pi/6)./x.^4; % number based distribution
    cryst_output.x_perc=x_perc;
    cryst_output.CSD_perc=CSD_perc/100; % volume percentage distribution 
    clear x

    % Create parameters object    
    p.min_length_discr=3e-4;   % grid spacing for deliquoring and drying
    p = carousel_parameters(p);     
    
    % time variables initialization
    process_time = 0;
    n_cycle = 0;
    batch_time = 0;    
   
    % other inizializations
    u_nominal=u;
    x=[];   
    measurements=[];
    y=[];
    
    %% Simulation
    while process_time <= total_duration 
        
        if batch_time < u.t_rot && n_cycle > 0   % simulate a step of process operation
            
           % duration of the next simulation step
           simulation_step = min(1, u.t_rot - batch_time); 
           
           % simulation
           [x,y,measurements]=carousel_simulator(batch_time,simulation_step,p,d,u,x,y,measurements,n_cycle);     
                     
           % update timers
           batch_time = batch_time + simulation_step; 
           process_time = process_time + simulation_step;
           
           % call online estimation routines
           x_estim = estimator_online(x_estim,u,y,control_flag);
           
           % call online control routines and save MVs profiles
           if ceil(batch_time/p.control_interval)== batch_time/p.control_interval
               [u,manipulated_vars] = controller_online(process_time,batch_time,...
               p.ports_working,u,u_nominal,measurements,manipulated_vars,x_estim,...
               n_cycle,control_flag);
           end
        else % rotation                      
           
           % call end of cycle estimation routines
           x_estim = estimator_cycle_switch(x_estim,u,y,control_flag);
            
           % calculate ethanol content in discharged cake
           y = final_composition(p,x,y,n_cycle); 
                                          
           % update timers
           batch_time = 0;
           n_cycle = n_cycle+1;
           p.ports_working=d.ports_working(n_cycle,:);
           
           % call end of cycle control routines and save MVs profiles
           [u,manipulated_vars] = controller_cycle_switch(process_time,batch_time,...
               p.ports_working,u,u_nominal,measurements,manipulated_vars,x_estim,...
               n_cycle,control_flag);
           
           % call disturbance function
           [cryst_output,d,p]=disturbances(process_time,cryst_output,p,d,u,n_cycle,disturbance_flag);

           % call switch cycle routine
           [x,y,measurements]=switch_cycle(process_time,cryst_output,p,d,u,x,y,measurements,n_cycle);

       end
    end

%% Prepare output object
if length(manipulated_vars.n_cycle_vector)>1
    d1.resistances=d.resistances(1:manipulated_vars.n_cycle_vector(end),:);
    d1.ports_working=d.ports_working(1:manipulated_vars.n_cycle_vector(end),:);
    d1.c_slurry=d.c_slurry(1:manipulated_vars.n_cycle_vector(end));
    d1.V_slurry=d.V_slurry(1:manipulated_vars.n_cycle_vector(end));
    d1.E=d.E(1:manipulated_vars.n_cycle_vector(end));
    d1.hM=d.hM(1:manipulated_vars.n_cycle_vector(end));
    d1.hT=d.hT(1:manipulated_vars.n_cycle_vector(end));
    d=d1;
else
    disp('No cakes discharged')
end

simulation_output.states=y;
simulation_output.measurements=measurements;
simulation_output.disturbances=d;
simulation_output.manipulated_vars=manipulated_vars;
simulation_output.estimated_states_parameters=x_estim;
simulation_output.feed.c_slurry_vector=cryst_output.conc_slurry;

simulation_output.settings.control_flag=control_flag;
simulation_output.settings.disturbance_flag=disturbance_flag;
simulation_output.settings.control_interval=control_interval;
simulation_output.settings.sampling_time=sampling_time;
simulation_output.settings.total_duration=total_duration;
simulation_output.settings.c_slurry_initial=cryst_output.conc_slurry;
simulation_output.settings.u_nominal=u_nominal;

if n_cycle>4
    throughput=d1.V_slurry.*simulation_output.manipulated_vars.V_slurry_vector.*...
        d1.c_slurry.*simulation_output.feed.c_slurry_vector;
    simulation_output.throughput=sum(throughput(1:length(simulation_output.states.final_composition)));
else
    simulation_output.throughput=0;
end
    simulation_output.manipulated_vars.t_rot_vector;