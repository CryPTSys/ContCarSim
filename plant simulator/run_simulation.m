function [process_time,p,x,y,controller_output]=run_simulation(p,u_ss,cryst_output,disturbance_flag,control_flag,cycles_number)  
    %% Initialization
    total_duration = u_ss.t_rot*cycles_number;
    load set_points    
    
    % load CSD
    load CSD
    cryst_output.x=x;
    cryst_output.CSD=CSD/100/(pi/6)./x.^4; % number based distribution
    cryst_output.x_perc=x_perc;
    cryst_output.CSD_perc=CSD_perc/100; % volume percentage distribution 
    clear x
    
    u_ss.flowrate_slurry=u_ss.V_slurry/u_ss.t_rot;% m3/s
    u=u_ss;     
    
    p.number_nodes_washing=50; % number of nodes washing (analytical solution)
    p.min_length_discr=4e-4;   % grid spacing for deliquoring and drying
    p = carousel_parameters(p);
    
    % time variables initialization
    process_time = 0;
    n_cycle = 0;
    cycle_time = 0;

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
        
        if cycle_time < u.t_rot % simulate a step of process operation
            
           % duration of the next simulation step
           simulation_step = min(p.control_interval, u.t_rot - cycle_time); 
           
           % simulation
           [x,y]=carousel_simulator(cycle_time,simulation_step,p,u,x,y,n_cycle);     
           
           % update timers
           cycle_time = cycle_time + simulation_step; 
           process_time = process_time + simulation_step;
           
           % call disturbance function
           [cryst_output,p]=disturbances(process_time,cryst_output,p,u_ss,disturbance_flag); 
           
           % call control routines and save MVs profiles
           u = controller(sp,cycle_time,u_ss,p,u,y,n_cycle,control_flag); % for DP and T_drying update
           controller_output.t_vector=[controller_output.t_vector process_time];
           controller_output.dP_vector=[controller_output.dP_vector u.dP_drying];
           controller_output.Tin_drying=[controller_output.Tin_drying u.Tinlet_drying];
           controller_output.rmv=[controller_output.rmv p.Rm];
           controller_output.t_rot_vector=[controller_output.t_rot_vector u.t_rot];
           
        else % rotation 
            
           % concatenate measurements of cycle that has just finished to measurements of previous cycles
           % for having continuous output from sensors
           y = continuous_outputs(process_time,p,x,y,n_cycle); 
           
           % call control routine calculating next cycle duration
           u = controller_cycle_duration(sp,cycle_time,u_ss,p,u,y,n_cycle,control_flag);

           % update timers
           cycle_time = 0;
           n_cycle = n_cycle +1;
           
           % call switch cycle routine
           [x,y]=switch_cycle(process_time,cryst_output,p,u,x,y,n_cycle);

        end 
    end
    
% concatenate measurements of last cycle to measurements of previous cycles  
% for having continuous output from sensors
y = continuous_outputs(process_time,p,x,y,n_cycle);



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
%
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