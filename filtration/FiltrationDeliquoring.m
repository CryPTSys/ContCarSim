%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtration-deliquoring model
% F. Destro, v1: March27, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc,close all

load CSD_daniel
CSD=CSD(end,:);
t=t(1):0.1:t(end);

%% Inputs
conc_MSMPR=11; % kg/m3   
flowrate_MSMPR=1.7e-7; % m3/s
t_rot=240; % s
dP=5e4; % Pa
CycleNo=1;%u(Nx+5);
t_process=t_rot;

%% Filtration/deliquoring parameters
Filter_d = 0.01;               % Filter_diameter [m]
Rm = 2.22e9;                   % Filter medium resistance [??m^-1];
rho_sol = 1400;                % Crystal density [kg/m^3]
rho_liq = 842;                 % Liquid density [kg/m^3]
visc = 1.4E-03;                % Fluid viscosity [Pas]
kappa = 1;                     % Dynamic shape factor, sphere= 1
Surf_t = 22.39e-3;             % Surface tension [N/m]   
delx=x(2)-x(1);

%% Calculations filtration (cake formation)
A = Filter_d.^2.*pi./4;                        % Filtration area [m^2]
E=Porosity_function(Filter_d,CSD,x);           % Ouchiyama model
eps_s= 1-E;                                    % Solid volume fraction in cake [m^3 cryst/m^3 cake]
m0=trapz(x,CSD); 
m1=trapz(x,CSD.*x'); 

V_slurry_initial=flowrate_MSMPR.*t_rot;        % Slurry volume fed at the beginning of the filtration batch [m^3]
%L_slurry_init=V_slurry_in./A;                 % Height of the slurry column at the beginning [m]

m_solid_initial=V_slurry_initial.*conc_MSMPR;    % Total solid mass filled into the filter [kg]
V_solid_initial=m_solid_initial./rho_sol;        % Total solid volume filled into the filter [m^3]
w_solid_initial=conc_MSMPR./(conc_MSMPR+rho_liq.*(1-conc_MSMPR./rho_sol));  % Solid mass fraction in the initial slurry [-]

mass_slurry_initial=m_solid_initial/w_solid_initial;            % Total cake mass filled into the filter [kg]

V_liq_initial=V_slurry_initial-V_solid_initial;  % Total liquid volume filled into the filter [m^3]

V_liquid_pores_filtration=V_solid_initial*E/(1-E); % at the end of filtration [m3]
V_filt_final=V_liq_initial-V_liquid_pores_filtration;   % volume of filtrate at the end of filtration [m3]

% c_brig= 1./(((1-w_solid_initial)./(w_solid_initial.*rho_liq))-...
%     ((E)./(eps_s.*rho_sol)));   % Mass of dry cake deposited per unit volume of filtrate [kg cryst/ m3 filtrate liquid]

c=(m_solid_initial)/V_filt_final; % Mass of dry cake deposited per unit volume of filtrate [kg sol/ m3 filtrate liquid]more compact way, but same results as Brigi's

% Cake resistance
alpha_CSD=@(dp) 180*(1-E)./(E^3*dp.^2*rho_sol);
% alpha= sum(alpha_CSD(x).*CSD')/m0*delx;
alpha=trapz(x,alpha_CSD(x).*CSD'/m0);

% Determine the filtration states as function of time                                       
t_filt_total=visc*alpha*c*V_filt_final^2/(2*A^2*dP)+visc*Rm*V_filt_final/(A*dP); 
t_filt=t(t<t_filt_total);  % filtration time
a_filt= visc.*alpha.*c./(2.*A.^2.*dP);  % darcy coeff 1
b_filt= visc.*Rm./(A.*dP);      % darcy coeff 2
V_filt= (-b_filt+sqrt(b_filt.^2+4.*a_filt.*t_filt))./(2.*a_filt);     % V(t) - filtrate volume [m^3]
m_dry_cake=c*V_filt;  % mass of deposited dry cake as function of time [kg]
m_wet_cake=m_dry_cake*(1+E/(1-E)/rho_sol*rho_liq);
Q=(alpha*visc*c*V_filt/(dP*A^2)+Rm*visc/(dP*A)).^-1;
H_cake=(m_dry_cake/rho_sol/(1-E))/A; % height of the deposited cake [m]
m_liq_pores=m_wet_cake-m_dry_cake; % mass of liquid in the pores during filtration [m3]
solvent_content_mass=m_liq_pores./(m_liq_pores+m_dry_cake);
solvent_content_vol_filt=m_liq_pores/rho_liq./(m_liq_pores/rho_liq+m_dry_cake/rho_sol);

%% Deliquoring  
t_deliq=linspace(0,t_rot-t_filt_total,round(t_rot-t_filt_total)*100);
L_cake=H_cake(end);                   % cake height at the end of filtration [m]
pb=(4.6.*(1-E).*Surf_t)./(E.*m1/m0);          % threshold vacuum [Pa]
k_av= 1./ (alpha.*rho_sol.*(1-E));      % cake permeability [?]
p_star=dP./pb;
N_cap=(E^3*(m1/m0)^2*(rho_liq*9.81.*L_cake+dP))/((1-E)^2*L_cake.*Surf_t);
S_inf= 0.155.*(1+0.031.*N_cap.^(-0.49));   % irreducible cake saturation (minimum saturation that can be achieved by displacement of the interstitial liquid by the applied vacuum

ThetaP=(t_deliq.*k_av.*dP)./(E.*visc.*(L_cake.^2).*(1-S_inf)); % Dimensionless deliquoring time [-]
n_switch=min(max(ThetaP-1.915,0)*10,1);           
SR=1./(1+1.46.*ThetaP.^0.48).*n_switch+1./(1+1.08.*ThetaP.^0.88).*(1-n_switch);

S=S_inf+SR.*(1-S_inf); 
V_deliq=(1-S)*V_liquid_pores_filtration;
V_liq_pores_deliq=S*V_liquid_pores_filtration;
m_liq_pores_deliq=V_liq_pores_deliq*rho_liq; % mass of liquid in the pores during filtration [m3]
m_liq_pores_eq=S_inf*V_liquid_pores_filtration*rho_liq; 
solvent_content_mass_deliq=m_liq_pores_deliq./(m_liq_pores_deliq+m_dry_cake(end));
solvent_content_vol_deliq=V_liq_pores_deliq./(V_liquid_pores_filtration+m_dry_cake(end)/rho_sol);
solvent_content_vol_eq=m_liq_pores_eq/rho_liq./(m_liq_pores_eq/rho_liq+m_dry_cake(end)/rho_sol);
scd=V_liq_pores_deliq./(V_liquid_pores_filtration+V_solid_initial);
%% Graphical output

% Solvent content during deliquoring
plot(t_deliq,solvent_content_vol_deliq,[t_deliq(1) t_deliq(end)],[solvent_content_vol_eq solvent_content_vol_eq],'linewidth',1.5)
xlabel('Deliquoring time [s]')
ylabel('Cake vol. solvent content [-]')
set(gca,'fontsize',16,'linewidth',1.3,'xlim',[t_deliq(1) 60],'xtick',0:10:60)
legend('Solvent content','Equilibrium moisture content')

% Plot CSD
figure
plot(x,CSD,'linewidth',1.3)
xlabel('Particle size [m]')
ylabel('f [#/m^4]')
set(gca,'fontsize',16,'linewidth',1.3) %,'xlim',[t_deliq(1) 60],'xtick',0:10:60)

% % Plot solvent content profile
% plot([t_filt t_deliq+t_filt(end)],[solvent_content_vol_filt solvent_content_vol_deliq],[0 t_filt(end)+t_deliq(end)],[solvent_content_vol_eq solvent_content_vol_eq],'linewidth',1.5)
% xlabel('Time [s]')
% ylabel('Cake vol. solvent content [-]')
% set(gca,'fontsize',16,'linewidth',1.3)%,'xlim',[t_deliq(1) 60],'xtick',0:10:60)
% legend('Solvent content','Equilibrium moisture content')

% Plot filtrate volume profile
figure
plot([t_filt t_deliq+t_filt(end)],[V_filt V_filt(end)+V_deliq],'linewidth',1.5)
hold on,plot([t_filt_total t_filt_total],[0 1e-4],'r','linewidth',1)
xlabel('Time [s]')
ylabel('Filtrate volume [m^3]')
set(gca,'fontsize',16,'linewidth',1.3,'ylim',[0 6e-5]) %,'xtick',0:10:60)

