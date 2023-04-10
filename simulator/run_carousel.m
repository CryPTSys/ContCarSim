%% Script for running the carousel simulator
% Francesco Destro, October2021

clc, clear, %close all
rng default

%% set simulation conditions    
control_mode= 0; % variable passed to controller.m function, useful to test multiple control strategies
                    % implemented control strategies:
                    % 0: open-loop
                    % 1: control strategy #3 in Destro et al. (2022); example of end-point controller on temperature (works with u_nominal.V_slurry=3e-6 and cryst_output.conc_slurry=250) 
disturbance_scenario = 0;  % 0: normal operating conditions; 1: nominal slurry concentration ramp; 2: cake resistance step
total_duration = 1800; % s

% set nominal operating variables   
u_nominal.t_cycle=30;           % nominal set-point of cycle duration (s) 
                              % u_nominal.t_cycle and u.t_cycle MUST ALWAYS BE INTEGERS
u_nominal.V_slurry=3e-6;     % nominal set-point of fed slurry volume (m3) 
                                  % Note: can't be larger than 10e-6; to process
                                  % larger slurry volumes, comment lines 15-17
                                  % of run_simulation.m
u_nominal.P_compr=10e4;            % nominal set-point of gauge pressure compressor (Pa)
u_nominal.Tinlet_drying=50+273.15;      % nominal set-point of drying gas temperature (K) 

% set nominal feed properties
cryst_output.conc_slurry=250;   % nominal slurry concentration (kg/m3) - 
                                % actual slurry concentration subject to Gaussian disturbances (+ ramp if disturbance_scenario==1)

%% set sampling interval and control time
control_interval = 1; % time step at which controller_online.m is called (s)
                        % MUST BE MULTIPLE Of 1 s                         
sampling_interval = .1; % sampling time for output measurements and states
                    % MUST BE SUBMULTIPLE OF 1 s

%% Set inter-cycle idle time and mesh cleaning idle time                    
inter_cycle_Dt = 0; % dead time at the end of every cycle (s); default = 0
mesh_clean_Dt  = 0; % dead time at mesh cleaning (s); default = 0
                    
%% run simulator
simulation_output=run_simulation(u_nominal,cryst_output,disturbance_scenario,...
    control_mode,total_duration,control_interval,sampling_interval,...
    inter_cycle_Dt,mesh_clean_Dt);

clearvars -except simulation_output 