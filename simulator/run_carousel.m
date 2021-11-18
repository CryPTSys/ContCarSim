%% Script for running the carousel simulator
% Francesco Destro, October2021

clc, clear, %close all
rng default

%% set simulation conditions    
control_mode= 1; % variable passed to controller.m function, useful to test multiple control strategies
                    % implemented control strategies:
                    % 0: open-loop
                    % 1: control strategy #3 in companion paper; example of end-point controller on temperature (works with u_nominal.V_slurry=3e-6) 
disturbance_scenario = 2;  % 0: normal operating conditions; 1: nominal slurry concentration ramp; 2: cake resistance step
total_duration = 1800; % s

%% set nominal manipulated variables   
u_nominal.t_rot=30;           % cycle duration (s) 
                              % u_nominal.t_rot and u.t_rot MUST ALWAYS BE INTEGERS
u_nominal.V_slurry=3e-6;      % fed slurry volume (m3) 
u_nominal.P_compr=10e4;            % gauge pressure compressor (Pa)
u_nominal.Tinlet_drying=50+273.15;      % drying gas temperature (K) 

%% set nominal feed properties
cryst_output.conc_slurry=250;   % initial nominal slurry concentration (kg/m3) - 
                                % subject to Gaussian disturbances (+ ramp if disturbance_flag==1)

%% set sampling interval and control time
control_interval = 1; % time step at which controller_online.m is called (s)
                        % MUST BE MULTIPLE Of 1 s                         
sampling_interval = .1; % sampling time for output measurements and states
                    % MUST BE SUBMULTIPLE OF 1 s

%% Set inter-cycle dead time and mesh cleaning duration                    
inter_cycle_Dt = 0; % dead time at the end of every cycle (s); default = 0
mesh_clean_Dt  = 0; % dead time at mesh cleaning (s); default = 0
                    
%% run simulator
simulation_output=run_simulation(u_nominal,cryst_output,disturbance_scenario,...
    control_mode,total_duration,control_interval,sampling_interval,...
    inter_cycle_Dt,mesh_clean_Dt);

clearvars -except simulation_output 