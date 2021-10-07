function [process_time,p,d,x,y,controller_output,x_estim]=run_simulation(p,u,...
    cryst_output,disturbance_flag,control_flag,total_duration,controller_output,x_estim)  
    
%% Initialization           
    % Initialize disturbances vector
    d.solids_dist=1; 
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
    cycle_time = 0;
    
    % measurements vector initialization
    y.pos1.cycle_1.t_filt=0;
    y.pos1.cycle_1.V_filt=0;
    y.cont_sign.pos1_4.bias=0;
    
    % other inizializations
    p.time_vector=[];
    p.Rm_array=[];
    p.n_cycle_cip=3;
    p.wg_end=[];
    x=[];   
    u.dP_drying=u.dP;     % pressure drop Station 4 [Pa]
        
    %% Simulation
    while process_time < total_duration 
        
        if cycle_time < u.t_rot && n_cycle >0   % simulate a step of process operation
            
           % duration of the next simulation step
           simulation_step = min(p.control_interval, u.t_rot - cycle_time); 
           
           % simulation
           [x,y]=carousel_simulator(cycle_time,simulation_step,p,d,u,x,y,n_cycle);     
                     
           % update timers
           cycle_time = cycle_time + simulation_step; 
           process_time = process_time + simulation_step;
           
           % call online estimation routines
           x_estim = estimator_online(x_estim,u,y,control_flag);
           
           % call online control routines and save MVs profiles
           [u,controller_output] = controller_online(process_time,u,y,controller_output,x_estim,n_cycle,control_flag);
           
        else % rotation                      
           
           % call end of cycle estimation routines
           x_estim = estimator_cycle_end(x_estim,u,y,control_flag);
            
           % call end of cycle control routines and save MVs profiles
           [u,controller_output] = controller_cycle_end(process_time,u,y,controller_output,x_estim,n_cycle,control_flag);
           
           % concatenate measurements of cycle that has just finished to measurements of previous cycles
           % for having continuous output from sensors
           y = continuous_outputs(process_time-cycle_time,p,x,y,n_cycle); 
          
           % update timers
           cycle_time = 0;
           n_cycle = n_cycle +1;
           p.ports_working=d.ports_working(n_cycle,:);
           
           % call disturbance function
           [cryst_output,d,p]=disturbances(process_time,cryst_output,p,d,u,n_cycle,disturbance_flag);

           % call switch cycle routine
           [x,y]=switch_cycle(process_time,cryst_output,p,d,u,x,y,n_cycle);
           
           % store previous cycle duration and slurry volume                      
           p.time_vector=[p.time_vector process_time];                  
                      
       end 
    end
    
% concatenate measurements of last cycle to measurements of previous cycles  
% for having continuous output from sensors
y = continuous_outputs(process_time,p,x,y,n_cycle);

d.resistances=d.resistances(1:controller_output.n_cycle_vector(end));