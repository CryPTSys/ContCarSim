%% Script for setting the operating conditions of the simulation
clc, clear
close all

% create object storing physical properties and operating conditions
p=carousel_parameters_class_control; 

% Set operating conditions    
disturbance_flag = 1;   % 0: NOC, 1: fouling, 2: increase of impurity in feed   
control_flag = 1;       % 0: open-loop, 1: PID control
cycles_number = 20;	

u_ss.t_rot=240; % s
u_ss.V_slurry=8e-6;
u_ss.W=1;               % Washing ratio
u_ss.dP=5e4;            % Pressure drop filtration, washing and deliquoring [Pa]
u_ss.dP_drying=u_ss.dP;     % Pressure drop drying [Pa]
u_ss.Tinlet_drying=70+273.15;       % Drying gas temperature [K]

cryst_output.conc_MSMPR=100;  % kg/m3    
cryst_output.liq_mass_fr_vect=[0.95 0 0.05]';  % 95% mother liquor, 5% impurity
cryst_output.T=298;
p.wash_solvent_mass_fr=[0 1 0]'; % mass fractions - components 1-3 

% Set sampling time and control time
p.control_interval = 10; % seconds
p.filtration_sampling_time = 1; % filtrate flowrate sampling time (positions 1-4)
p.drying_sampling_time = 1; % gas temperature and composition sampling time (position 4)

% run simulator
[t,p,x,y,controller_output]=run_simulation(p,u_ss,cryst_output,disturbance_flag,...
    control_flag,cycles_number);

% %% set-point storage
% gas temperature
sp.Tg_pos4_ref=y.pos4.cycle_1.Tg;
sp.t_ref_Tg=y.pos4.cycle_1.t_drying;
sp.t_rot_ref=u_ss.t_rot;
% 
% % filtrate volume
% [~,starting]=min(abs(y.cont_sign.pos1_4.t-u_ss.t_rot*5));
% [~,ending]=min(abs(y.cont_sign.pos1_4.t-u_ss.t_rot*6));
% sp.filtrate_pos14_ref=y.cont_sign.pos1_4.V(starting+1:ending);
% sp.t_ref_filt=y.cont_sign.pos1_4.t(starting+1:ending)-u_ss.t_rot*5;
% sp.t_rot_ref=u_ss.t_rot; %s
% sp.V_slurry=u_ss.V_slurry; % m2
% save('set_points','sp')

%% graphical output
figure(1)
plot(p.time_vector/60, p.Rm_vector, 'k','linewidth',1)
xlabel('Time [min]')
ylabel('Filter mesh resistance [1/m]')
set(gca,'fontsize',18,'linewidth',1)
box on

figure(2)
box on
hold on
plot(linspace(0,t,length(controller_output.dP_vector))/60,controller_output.t_rot_vector/60,'linewidth',1.5)
set(gca,'fontsize',18,'linewidth',1)
xlabel('Time [min]')
ylabel('Cycle duration [min]')
% legend('slurry conc = 50 kg/m3 - nominal','slurry conc = 90 kg/m3','slurry conc = 10 kg/m3 ')

figure(3)
box on
semilogy(max(y.final_composition(1:end-1,:)'),'linewidth',1)
hold on
set(gca,'fontsize',18,'linewidth',1)
xlabel('Cycle #')
ylabel('Max impurity vol. fract. [-]')
% legend('slurry conc = 50 kg/m3 - nominal','slurry conc = 90 kg/m3','slurry conc = 10 kg/m3 ')

figure(4)
plot(y.pos4.cycle_1.t_drying,y.pos4.cycle_1.Tg),hold on
plot(y.pos4.cycle_5.t_drying,y.pos4.cycle_5.Tg)
plot(y.pos4.cycle_10.t_drying,y.pos4.cycle_10.Tg)
plot(y.pos4.cycle_15.t_drying,y.pos4.cycle_15.Tg)
