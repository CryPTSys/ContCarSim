%% Script for running the carousel simulator
% Francesco Destro, October2021

clc, clear, %close all
rng default

%% set simulation conditions    
control_mode= 2; % variable passed to controller.m function, useful to test multiple control strategies
                    % implemented control strategies:
                    % 0: open-loop
                    % 1: control strategy #3 in companion paper; example of end-point controller on temperature (works with u_nominal.V_slurry=3e-6) 
disturbance_scenario = 0;  % 0: normal operating conditions; 1: nominal slurry concentration ramp; 2: cake resistance step
total_duration = 3600; % s

%% set nominal operating variables   
u_nominal.t_cycle=30;           % cycle duration (s) 
                              % u_nominal.t_cycle and u.t_cycle MUST ALWAYS BE INTEGERS
u_nominal.V_slurry=6e-6;      % fed slurry volume (m3) 
u_nominal.P_compr=10e4;            % gauge pressure compressor (Pa)
u_nominal.Tinlet_drying=50+273.15;      % drying gas temperature (K) 

%% set nominal feed properties
cryst_output.conc_slurry=250;   % initial nominal slurry concentration (kg/m3) - 
                                % subject to Gaussian disturbances (+ ramp if disturbance_flag==1)

%% set sampling interval and control time
control_interval = 1; % time step at which controller_online.m is called (s)
                        % MUST BE MULTIPLE Of 1 s                         
sampling_interval = .2; % sampling time for output measurements and states
                    % MUST BE SUBMULTIPLE OF 1 s

%% Set inter-cycle idle time and mesh cleaning idle time                    
inter_cycle_Dt = 0; % dead time at the end of every cycle (s); default = 0
mesh_clean_Dt  = 0; % dead time at mesh cleaning (s); default = 0
                    
%% run simulator
simulation_output=run_simulation(u_nominal,cryst_output,disturbance_scenario,...
    control_mode,total_duration,control_interval,sampling_interval,...
    inter_cycle_Dt,mesh_clean_Dt);

clearvars -except simulation_output 