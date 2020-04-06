%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtration-deliquoring model
% F. Destro, v1: March27, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
% clc,close all

%% Load data from file, in the future these data will be provided by cryst model
load CSD_daniel
CSD=CSD(end,:);
cryst_output.x=x;
cryst_output.CSD=CSD;
cryst_output.conc_MSMPR=11; % kg/m3   
cryst_output.flowrate_MSMPR=1.7e-7; % m3/s
clear x, clear CSD

%% Inputs
p.t_rot=200; % s
p.dP=5e4; % Pa
t=0:0.01:p.t_rot;
% CycleNo=1;%u(Nx+5);
% t_process=3600;
% t=t(1):0.1:t_process);

%% Filtration/deliquoring parameters
p.Filter_d = 0.01;               % p.Filter_diameter [m]
p.Rm = 2.22e9;                   % Filter medium resistance [??m^-1];
p.rho_sol = 1400;                % Crystal density [kg/m^3]
p.rho_liq = 842;                 % Liquid density [kg/m^3]
p.visc = 1.4E-03;                % Fluid p.viscosity [Pas]
p.kappa = 1;                     % Dynamic shape factor, sphere= 1
p.surf_t = 22.39e-3;             % Surface tension [N/m]   

%% Simulation section

% Filtration 
filt_output=model_filtration(t,cryst_output,p);

% Deliquoring
% Select deliquoring length: up to t_rot or further?
p.t_deliq_final=p.t_rot-filt_output.t_filt_total; 
% % Solve with design charts and eq. conc. solvent depending on Dmean
deliq_output=model_deliquoring_design_charts(filt_output,p);
% Solve with design charts and eq. conc. solvent depending on CSD
deliq_output2=model_deliquoring_design_charts_eqCSD(cryst_output,filt_output,p);
% Solve integrating the PDEs

%% Graphical output

t_deliq=deliq_output.t_deliq;
if length(t_deliq)<1
    t_deliq=0;
    V_deliq=0;
    solvent_content_vol_deliq=filt_output.solvent_content_vol_filt(end);
    solvent_content_vol_eq=deliq_output.solvent_content_vol_eq;
else
    V_deliq=deliq_output.V_deliq;
    solvent_content_vol_deliq=deliq_output.solvent_content_vol_deliq;
    solvent_content_vol_eq=deliq_output.solvent_content_vol_eq;
    S=deliq_output.S;
    S_inf=deliq_output.S_inf;
end
CSD=cryst_output.CSD;
x=cryst_output.x;
t_filt=filt_output.t_filt;
V_filt=filt_output.V_filt;
solvent_content_vol_filt=filt_output.solvent_content_vol_filt;
t_filt_total=filt_output.t_filt_total;

% % Solvent content during deliquoring
% plot(t_deliq,solvent_content_vol_deliq,[t_deliq(1) t_deliq(end)],[solvent_content_vol_eq solvent_content_vol_eq],'linewidth',1.5)
% xlabel('Deliquoring time [s]')
% ylabel('Cake vol. solvent content [-]')
% set(gca,'fontsize',16,'linewidth',1.3,'xlim',[t_deliq(1) t_deliq(end)],'xtick',round(linspace(0,t_deliq(end),6)))
% legend('Solvent content','Equilibrium moisture content')

% Plot CSD
% figure
% plot(x,CSD,'linewidth',1.3)
% xlabel('Particle size [m]')
% ylabel('f [#/m^4]')
% set(gca,'fontsize',16,'linewidth',1.3) %,'xlim',[t_deliq(1) 60],'xtick',0:10:60)

% Plot solvent content profile during filtration and deliquoring
figure
plot([t_filt t_deliq+t_filt(end)],[solvent_content_vol_filt solvent_content_vol_deliq],...
    [0 t_filt(end)+t_deliq(end)],[solvent_content_vol_eq solvent_content_vol_eq],'linewidth',1.5);
lim=get(gca);
lim=lim.YLim;
hold on,plot([t_filt_total t_filt_total],[0 lim(2)],'r','linewidth',1)
xlabel('Time [s]')
ylabel('Cake vol. solvent content [-]')
set(gca,'fontsize',16,'linewidth',1.3)%,'xlim',[t_deliq(1) 60],'xtick',0:10:60)
legend('Solvent content','Equilibrium moisture content')
% 
% % Plots saturation profile
% figure
% plot(t_deliq,S,[0 t_deliq(end)],[S_inf S_inf],'linewidth',1.5)
% xlabel('Time [s]')
% ylabel('Cake saturation [-]')
% set(gca,'fontsize',16,'linewidth',1.3)%,'xlim',[t_deliq(1) 60],'xtick',0:10:60)
% legend('Saturation','Equilibrium saturation')
% 
% 
% % Plot filtrate profile
% figure
% plot([t_filt t_deliq+t_filt(end)],[V_filt V_filt(end)+V_deliq],'linewidth',1.5)
% hold on,plot([t_filt_total t_filt_total],[0 1e-4],'r','linewidth',1)
% xlabel('Time [s]')
% ylabel('Filtrate volume [m^3]')
% set(gca,'fontsize',16,'linewidth',1.3,'ylim',[0 6e-5]) %,'xtick',0:10:60)

%% Warning
if length(deliq_output.ThetaP)>1
    if deliq_output.ThetaP(end)>204
        disp('Attention! The solution has been obtained through extrapolation')
    end
end