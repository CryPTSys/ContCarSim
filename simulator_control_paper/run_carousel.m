%% Script for running the carousel simulator
clc, clear, close all

rng default

%% set simulation conditions    
control_flag = 0; % variable passed to controller.m function, useful to test multiple control strategies
disturbance_flag = 0;   % 0: only fouling, 1: fouling + slurry conc. ramp
total_duration = 600;  % s

%% set nominal manipulated variables   
u_ss.t_rot=75;           % cycle duration (s) 
u_ss.V_slurry=6e-6;      % fed slurry volume (m3)
u_ss.dP=10e4;            % pressure drop stations 1-4 (Pa)
u_ss.Tinlet_drying=50+273.15;       % drying gas temperature (K)

%% set nominal feed properties
cryst_output.conc_MSMPR=250;   % (kg/m3)
cryst_output.T=295.25;         % (K)

%% set sampling time and control time
p=carousel_parameters_class_control; % do not edit - create object storing physical properties and operating conditions
p.control_interval = 1; % time step at which control routine is called (s)
                        % SET as multiple of p.filtration_sampling_time and
                        % p.drying_sampling_time
p.filtration_sampling_time = .2; % Stations 1-4 sensors sampling time
p.drying_sampling_time = .2; % Station 5 sensors sampling time

%% initialize vector containing MVs profiles updated by controller.m (other fields can be added)
% updated online
controller_output.t_vector=[]; % process time vector
controller_output.dP_vector=[];
controller_output.Tin_drying_vector=[];
% updated at the end of every cycle
controller_output.n_cycle_vector=[];
controller_output.t_rot_vector=[];  % rotation time vector
controller_output.V_slurry_vector=[]; % fed slurry vector

%% initialize estimated states/parameters vector (add fields)
x_estim=[]; % add fields if you set up a state estimator in controller.m

%% run simulator
[t,p,d,x,y,controller_output,x_estim]=run_simulation(p,u_ss,cryst_output,disturbance_flag,...
    control_flag,total_duration,controller_output,x_estim);