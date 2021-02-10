function [x,y,controller_output]=run_simulation()
    %% create object storing physical properties and operating conditions
    p=carousel_parameters_class_control; 
    
    %% Set operating conditions
    p.wash_solvent_mass_fr=[0 1 0]'; % mass fractions - components 1-3 
    disturbance_flag = 0; % 0: NOC, 1: fouling, 2: increase of impurity in feed   
    cycles_number = 10;	
    
    u_ss.t_rot=120; % s
    u_ss.V_slurry=5e-6;
    u_ss.W=1; % Washing ratio
    u_ss.dP=5e4; % Pressure drop filtration, washing and deliquoring [Pa]
    u_ss.dP_drying=5e4; % Pressure drop drying [Pa]
    u_ss.Tinlet_drying=70+273.15; % Drying gas temperature [K]
    
    cryst_output.conc_MSMPR=100;  % kg/m3   
    cryst_output.liq_mass_fr_vect=[0.95 0 0.05]';  % 99% mother liquor, 1% impurity
    cryst_output.T=298;
    
    % Set sampling time and control time
    p.control_interval = 10; % seconds
    p.filtration_sampling_time = 1; % filtrate flowrate sampling time (positions 1-4)
    p.drying_sampling_time = .1; % gas temperature and composition sampling time (position 4)
    
    %% Initialization
    total_duration = u_ss.t_rot*cycles_number;
    load set_points    
    
    load CSD_old
    CSD=CSD(end,:);
    CSD=CSD/sum(CSD);
    CSD=CSD(end,:);
    cryst_output.x=x;
    cryst_output.CSD=CSD';
    clear x, clear CSD       
    
    u_ss.flowrate_slurry=u_ss.V_slurry/u_ss.t_rot;%.025e-6;%1.7e-7; % m3/s
    u=u_ss;     
    
    p.number_nodes_washing=50; % number of nodes washing (analytical solution)
    p.min_length_discr=2e-4;   % grid spacing for deliquoring and drying
    p = carousel_parameters(p);
    
    % time variables initialization
    process_time = 0;
    n_rotation = 0;
    rotation_time = 0;

    % initial states 
    x.pos0.charge_cell_volume = 0;
    x.slurry_tank_volume = 20e-6;   
    
    % measurements vector initialization
    y.pos1.cycle_1.t_filt=0;
    y.pos1.cycle_1.V_filt=0;
    
    % initialize vector containing MVs profiles
    controller_output.dP_vector=[];
    controller_output.rmv=[];
    controller_output.t_vector=[];
    controller_output.Tin_drying=[];
    controller_output.t_rot_vector=[];
    
    u.error_t_rot(1:5)=0;
    p.time_vector=[];
    p.Rm_vector=[];
    
    %% Simulation
    while process_time < total_duration
        if rotation_time < u.t_rot % simulate a step of process operation
           % duration of the next simulation step
           step = min(p.control_interval, u.t_rot - rotation_time); 
           
           % simulation
           [x,y]=carousel_simulator(rotation_time,step,cryst_output,p,u,x,y,n_rotation);     
           
           % update timers
           rotation_time = rotation_time + step; 
           process_time = process_time + step;
           
           % call disturbance function
           [cryst_output,p,u]=disturbances(process_time,cryst_output,p,u,disturbance_flag); 
           
           % call control routines and save MVs profiles
           %u = controller(sp,rotation_time,u_ss,p,u,y,n_rotation); % for DP and T_drying update
           controller_output.t_vector=[controller_output.t_vector process_time];
           controller_output.dP_vector=[controller_output.dP_vector u.dP_drying];
           controller_output.Tin_drying=[controller_output.Tin_drying u.Tinlet_drying];
           controller_output.rmv=[controller_output.rmv p.Rm];
           controller_output.t_rot_vector=[controller_output.t_rot_vector u.t_rot];
           
        else % rotation 
           % concatenate measurements of new cycle to measurements of previous cycles
           y = continuous_outputs(process_time,p,x,y,n_rotation); 
           
           % call control routine calculating next cycle duration
%            u = controller_slow(sp,rotation_time,u_ss,p,u,y,n_rotation);

           % update timers
           rotation_time = 0;
           n_rotation = n_rotation +1;
           
           % call switch phase routine
           [x,y]=switch_phase(process_time,cryst_output,p,u,x,y,n_rotation);

        end 
    end
    
% concatenate measurements of last cycle to measurements of previous cycles  
y = continuous_outputs(process_time,p,x,y,n_rotation);


%% Graphical output
% figure(1) 
% box on
% hold on
% plot(y.pos4.cycle_2.t_drying,y.pos4.cycle_2.Tg,'linewidth',1.5)
% set(gca,'fontsize',18,'linewidth',1,'xtick',0:30:180,'xlim',[0 180])
% xlabel('Time [s]')
% ylabel('Drying gas T_{out} [K]')
% legend('slurry conc = 50 kg/m3 - nominal','slurry conc = 90 kg/m3','slurry conc = 10 kg/m3 ')

% figure(2)
% box on
% hold on
% plot(linspace(0,total_duration,length(controller_output.dP_vector)),controller_output.dP_vector,'linewidth',1.5)
% set(gca,'fontsize',18,'linewidth',1)
% xlabel('Time [s]')
% ylabel('Drying pressure drop [Pa]')
% legend('slurry conc = 50 kg/m3 - nominal','slurry conc = 90 kg/m3','slurry conc = 10 kg/m3 ')
% 
% figure(3)
% box on
% hold on
% plot(linspace(0,total_duration,length(controller_output.dP_vector)),controller_output.Tin_drying,'linewidth',1.5)
% set(gca,'fontsize',18,'linewidth',1)
% xlabel('Time [s]')
% ylabel('Drying inlet T [K]')
% legend('slurry conc = 50 kg/m3 - nominal','slurry conc = 90 kg/m3','slurry conc = 10 kg/m3 ')
figure(1)
plot(p.time_vector/60, p.Rm_vector, 'k','linewidth',1)
xlabel('Time [min]')
ylabel('Filter mesh resistance [1/m]')
set(gca,'fontsize',18,'linewidth',1)

figure(4)
box on
hold on
plot(linspace(0,total_duration,length(controller_output.dP_vector))/60,controller_output.t_rot_vector/60,'linewidth',1.5)
set(gca,'fontsize',18,'linewidth',1)
xlabel('Time [min]')
ylabel('Cycle duration [min]')
% legend('slurry conc = 50 kg/m3 - nominal','slurry conc = 90 kg/m3','slurry conc = 10 kg/m3 ')

figure(5)
box on
semilogy(max(y.final_composition(1:end-1,:)'),'linewidth',1)
hold on
set(gca,'fontsize',18,'linewidth',1)
xlabel('Cycle #')
ylabel('Max impurity vol. fract. [-]')
% legend('slurry conc = 50 kg/m3 - nominal','slurry conc = 90 kg/m3','slurry conc = 10 kg/m3 ')

%% set-point calculation
% % gas temperature
% sp.Tg_pos4_ref=y.pos4.cycle_1.Tg;
% sp.t_ref_Tg=0:p.drying_time_step:u.t_rot;
% sp.t_rot_ref=u.t_rot;
% save('set_points','sp')
% 
% % filtrate volume
% [~,starting]=min(abs(y.cont_sign.pos1_4.t-u.t_rot*5));
% [~,ending]=min(abs(y.cont_sign.pos1_4.t-u.t_rot*6));
% sp.filtrate_pos14_ref=y.cont_sign.pos1_4.V(starting+1:ending);
% sp.t_ref_filt=y.cont_sign.pos1_4.t(starting+1:ending)-u.t_rot*5;
% sp.t_rot_ref=u_ss.t_rot; %s
% sp.V_slurry=u_ss.V_slurry; % m2
% save('set_points','sp')

%     plot(t_vector,rmv)
%     set(gca,'fontsize',16,'linewidth',1)
%     xlabel('Time [s]')
%     ylabel('Filter resistance [1/m]')
%     
%     figure
%     plot(t_vector,dP_vector)
%     set(gca,'fontsize',16,'linewidth',1)
%     xlabel('Time [s]')
%     ylabel('Pressure drop dryer [-]')
% plot(y.pos4.cycle_1.t_drying,y.pos4.cycle_1.Tg)
% hold on
% plot(y.pos4.cycle_1.t_drying,y.pos4.cycle_14.Tg)
% set(gca,'fontsize',16,'linewidth',1)
% xlabel('Drying time [s]')
% ylabel('Dryer outlet temperature [-]')
% legend('NOC - sp','Batch starting @t=195 s, OL',...
%    'Batch starting @t=195 s, CL')
% 
% plot(y.final_composition(:,2),'o')
% set(gca,'fontsize',16,'linewidth',1)
% xlabel('Batch #')
% ylabel([{'Final vol. content of'},{'most concentrated impurity'} ])