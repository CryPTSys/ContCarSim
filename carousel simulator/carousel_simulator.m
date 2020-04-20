%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtration-deliquoring model
% F. Destro, v1: March27, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc,close all

%% Load data from file, in the future these data will be provided by cryst model
load CSD_daniel
CSD=CSD(end,:);
cryst_output.x=x;
cryst_output.CSD=CSD;
cryst_output.conc_MSMPR=11; % kg/m3   
cryst_output.flowrate_MSMPR=1.7e-7; % m3/s
clear x, clear CSD

%% Inputs
p.t_rot=240; % s
p.dP=5e4; % Pa
p.W=5; % washing ratio
t=0:0.01:p.t_rot;
% CycleNo=1;%u(Nx+5);
% t_process=3600;
% t=t(1):0.1:t_process);

%% Filtration/deliquoring/washing parameters
p.Filter_d = 0.01;               % p.Filter_diameter [m]
p.Rm = 2.22e9;                   % Filter medium resistance [??m^-1];
p.rho_sol = 1400;                % Crystal density [kg/m^3]
p.rho_liq = 842;                 % Liquid density [kg/m^3]
p.visc = 1.4E-03;                % Fluid viscosity [Pa s]
p.viscG = 1.8e-5;                % Air viscosity at room temperature and 1 atm [Pa s]
p.surf_t = 22.39e-3;             % Surface tension [N/m]  
p.lambda = 5;                    % Pore size index for deliquoring [-]
% p.grid_deliq = 1e-4;             % Grid spacing for deliquoring [m]
% p.grid_washing = p.grid_deliq;
p.number_nodes_deliq=50;
p.number_nodes_washing=50;
p.wash_time_step=1e-4;
p.Di = 1.28e-9;                  % diffusivity of ethanol in water [m2/s]

%% Calculation of cake properties 
p.A = p.Filter_d.^2.*pi./4;                        % Filtration area [m^2]
p.E=Porosity_function(p.Filter_d,cryst_output.CSD,cryst_output.x);           % Ouchiyama model
p.m0=trapz(cryst_output.x,cryst_output.CSD); 
p.m1=trapz(cryst_output.x,cryst_output.CSD.*cryst_output.x'); 
p.alpha_CSD=@(dp) 180*(1-p.E)./(p.E^3*dp.^2*p.rho_sol);
p.alpha=trapz(cryst_output.x,p.alpha_CSD(cryst_output.x).*cryst_output.CSD'/p.m0);
p.k= 1/(p.alpha*p.rho_sol*(1-p.E));      % cake permeability [?]
%% Simulation section

% Filtration 
[filt_output,p]=model_filtration(t,cryst_output,p);

% Deliquoring
% Select deliquoring length: up to t_rot or further?
p.t_deliq_final=p.t_rot-filt_output.t_filt_total; 
% % Solve with design charts  eq. conc. of solvent and treshold pressure depend on Dmean
% deliq_output=model_deliquoring_design_charts(filt_output,p);
% Solve with design charts - eq. conc. of solvent and treshold pressure depend on CSD
deliq_output_charts=model_deliquoring_design_charts(p);
% Solve integrating the PDEs - eq. conc. of solvent and treshold pressure depend on CSD
deliq_output_pde=model_deliquoring_pde_adim(p);

% Washing
% washing_output_pre_deliq_approx=model_washing_multiple_saturat_approx(deliq_output_pde,p);
washing_output=model_washing(deliq_output_pde,p);

%% Graphical output
CSD=cryst_output.CSD;
x=cryst_output.x;
t_filt=filt_output.t_filt;
V_filt=filt_output.V_filt;
solvent_content_vol_filt=filt_output.solvent_content_vol_filt;
t_filt_total=filt_output.t_filt_total;

t_deliq=deliq_output_charts.t_deliq;
if length(t_deliq)<1
    t_deliq=0;
%     V_deliq=0;
    solvent_content_vol_deliq_charts=filt_output.solvent_content_vol_filt(end);
    
else
%     V_deliq=deliq_output_charts.V_deliq;
    solvent_content_vol_deliq_charts=deliq_output_charts.solvent_content_vol_deliq;
    solvent_content_vol_deliq_pde=deliq_output_pde.solvent_content_vol_deliq;
    S=deliq_output_charts.S;
    
end

solvent_content_vol_eq=p.S_inf*p.E;

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
plot([t_filt t_deliq+t_filt(end)],[solvent_content_vol_filt solvent_content_vol_deliq_charts],...
    [t_filt t_deliq+t_filt(end)],[solvent_content_vol_filt solvent_content_vol_deliq_pde],...
    [0 t_filt(end)+t_deliq(end)],[solvent_content_vol_eq solvent_content_vol_eq])

xlabel('Time [s]')
ylabel('Cake mean vol. solvent content [-]')
set(gca,'fontsize',16,'linewidth',1.3,'xtick',0:20:500)%'xlim',[t_deliq(1) 60],'xtick',0:10:60)
axis([0 t_filt(end)+t_deliq(end) 0 solvent_content_vol_filt(end)*1.2] )
lim=get(gca);
lim=lim.YLim;
hold on,plot([t_filt_total t_filt_total],[0 lim(2)],'r','linewidth',.5)
legend('Solvent content - design charts','Solvent content - PDE','Equilibrium moisture content','Beginning deliquoring step')

figure
plot(deliq_output_pde.nodes_deliq,deliq_output_pde.S_final*p.E)
xlabel('Cake axial coordinate')
ylabel('Cake vol. solvent content profile end of deliquoring [-]')
set(gca,'fontsize',16,'linewidth',1.3)

% Plot solvent concentration after washing
figure
plot(deliq_output_pde.nodes_deliq,washing_output.xv_mother_liquor)
xlabel('Cake axial coordinate')
ylabel('Cake vol. solvent content profile end of washing[-]')
set(gca,'fontsize',16,'linewidth',1.3)

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
if length(deliq_output_charts.ThetaP)>1
    if deliq_output_charts.ThetaP(end)>204
        disp('Attention! The solution obtained through design charts is extrapolated')
    end
end
