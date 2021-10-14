%% Script for running the carousel simulator
% Francesco Destro, October2021

clc, clear, %close all
rng default

%% set simulation conditions    
control_flag = 0;   % variable passed to controller.m function, useful to test multiple control strategies
                    % implemented control strategies: 0: open-loop; 1: end-loop controller on temperature (works with u_nominal.V_slurry=6e-6)  
disturbance_flag = 0;  % 0: no disturbance; 1: nominal slurry concentration ramp
total_duration = 200;  % s

%% set nominal manipulated variables   
u_nominal.t_rot=1;            % cycle duration (s) 
u_nominal.V_slurry=5e-6;      % fed slurry volume (m3) - Range: [0.5e-6, 10e-6]
u_nominal.dP=10e4;            % pressure drop stations 1-4 (Pa)
u_nominal.Tinlet_drying=50+273.15;       % drying gas temperature (K)

%% set nominal feed properties
cryst_output.conc_slurry=250;   % initial nominal slurry concentration (kg/m3) - 
                                % subject to Gaussian disturbances (+ ramp if disturbance_flag==1)
cryst_output.T=295.25;          % (K)

%% set sampling time and control time
control_interval = 1; % time step at which control routine is called (s)
                        % SET as multiple of 1 s                         
sampling_time = .25; % Stations 1-5 sensors sampling time

%% initialize vector containing MVs profiles updated by controller.m (other fields can be added)
% updated online
manipulated_vars.t_vector=[]; % process time vector
manipulated_vars.dP_vector=[];
manipulated_vars.Tin_drying_vector=[];
% updated at the end of every cycle
manipulated_vars.n_cycle_vector=[];
manipulated_vars.t_rot_vector=[];  % rotation time vector
manipulated_vars.V_slurry_vector=[]; % fed slurry vector

%% initialize estimated states/parameters vector (add fields)
x_estim=[]; % add fields if you set up a state/parameter estimator in estimator_online.m and/or estimator_cycle_switch.m

%% run simulator
simulation_output=run_simulation(u_nominal,cryst_output,disturbance_flag,...
    control_flag,total_duration,manipulated_vars,x_estim,control_interval,sampling_time);

clearvars -except simulation_output 