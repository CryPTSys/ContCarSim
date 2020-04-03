%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulator for the convective dryer model %
% v1: F. Destro, March2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc

%% Parameters
% simulation inputs
p.n_nodes = 30;
p.tf = 5*8000; % process duration [s]
graph_interval=5000;
step_mom = 1000; % time step for recalculating the pressure profile
p.Z = 0.05; % dryer length [m]
p.D_dryer = 0.01; % dryer diameter [m]
p.eps_s = 0.6; % 1 - porosity [m3_s/m3]
p.epsL0 = 0.10; % volumetric fraction

% adjustable parameters
p.Tg_inlet = 70+273.15; % BC - [K] paracetamol melts at 440 K
p.Vg_N = 17e-6; % flowrate of drying nitrogen in [Nm3/s] - ref: 10NmL/s
p.d_p = 1e-4; % particles mean diameter [m]
p.a_V = 3500; %6*(1-p.eps_s)/p.d_p/5; % specific area [m2/m3]
p.eps_l_crit = 0.05; % critical solvent content [m3_l/m3]
p.h_M = 5e-10; % solvent mass transfer coefficient [1/s?] - model it!
        % from literature: as long as the driving force is positive, the
        % mass transfer occurs rapidly

% fixed parameters
p.eps_l_eq = 0; % equilibrium solvent content [m3_l/m3]
p.f=@(x) min((x(:,2)-p.eps_l_eq)/(p.eps_l_crit-p.eps_l_eq),1); % drying rate efficiency function [adim]
p.rho_s = 1400; % solid density [kg_s/m3] - aspirin
% p.rho_s = 1300; % solid density [kg_s/m3] - paracetamol
p.A_dryer = pi/4*p.D_dryer^2; % dryer cross-section [m2]
% p.cp_s = 190/150*1e3; % solid specific heat [J/(kg K)] paracetamol
p.cp_s = 2200; % aspirin % solid specific heat [J/(kg K)]
%p.h_T = 30; % heat transfer coefficient [1/s?] - model it!
p.A1 = 8.10765; % Antoine constant #1 % 0 to 60 C, Felder Rousseau
p.B1 = 1750.286; % Antoine constant #2
p.C1 = 235; % Antoine constant #3
p.A2 = 7.96681; % 60 to 150 C
p.B2 = 1668.21; 
p.C2 = 228; 
p.k_N2 = 0.028; % air conductivity [W/(m K)]
p.Dwat_air = 0.4e-4; % diffusivity of water in air [m2/s]
p.rho_l = 1000; % solvent liquid density [kg_l/m3]
p.lambda = 2260*1e3; % solvent latent heat of vaporization J/kg] - model it
p.mu = 1.76e-5; % gas phase viscosity [Pa s] - model it
p.MW_N2 = 28*1e-3; % air molecular weight [kg/mol]
p.MW_solv = 18*1e-3; % solvent molecular weight [kg/mol]
p.cp_solventG = [33.46e-3;0.688e-5;0.7604e-8;-3.593e-12]/p.MW_solv*1e3; % solvent (steam) specific heat in the gas phase [J/(kg K)] 
p.cp_N2 = [29e-3;0.2199e-5;0.5723e-8;-2.871e-12]/p.MW_N2*1e3; % nitrogen specific heat [J/(kg K)] 
p.cp_l = 4200; % liquid specific heat [J/(kg K)] 
% p.Pg_inlet = 101325*10; % BC - [Pa] - not used, the inlet pressure is calculated with the Momentum Balance
p.Pg_outlet = 101325; % BC - [Pa]

% methods
p.switch_antoine=@(x) 1-min(max(x(:,3)-60-273.15,0),1);
p.coeff= @(x) [p.A1,p.B1,p.C1].*p.switch_antoine(x)+[p.A2,p.B2,p.C2].*(1-p.switch_antoine(x)); 
p.cp = @(x,coeff) (coeff(1)+coeff(2)*x(:,3)+coeff(3)*x(:,3).^2+coeff(4)*x(:,3).^3); % [ J/kg K]
p.eps_g = @(x) 1-p.eps_s-x(:,2);
p.MWgas = @(x) (1-x(:,1))*p.MW_N2+x(:,1)*p.MW_solv;
p.Vg = @(x) p.Vg_N * (101325./x(:,5)).*(x(:,3)/273.15); % flowrate of drying nitrogen in [m3/s]
p.ug = @(x) p.Vg(x)/p.A_dryer./p.eps_g(x); % gas velocity [m/s]
p.rho_g=@(x) x(:,5)./(8.314*x(:,3)).*p.MWgas(x); % gas phase density [kg/m3]
p.cp_gas = @(x) p.cp(x,p.cp_N2).*(1-x(:,1))+p.cp(x,p.cp_solventG).*x(:,1);
p.Psat_solv = @(x,coeff) 10.^(coeff(:,1)-coeff(:,2)./(coeff(:,3)+(x(:,3)-273.15)))*133.322; % saturation pressure [Pa] 

p.Sc= @(x) p.mu/(p.rho_g(x)*p.Dwat_air);
p.Pr= @(x) p.cp_gas(x)*p.mu/p.k_N2;
p.Re = @(x) p.rho_g(x).*p.ug(x).*p.eps_g(x)*p.d_p/p.mu;
p.Re_bird = @(x) p.Re(x)/p.eps_s;
p.Nu = @(x) (2+1.3*p.Pr(x).^0.15+0.66*p.Pr(x).^0.31.*p.Re(x).^0.50);%.5; % Re is very low, so we can assume pure conduction - 0.79*(2.06/p.eps_g(x)*p.Re(x)^0.425*p.Pr(x)^(1/3));
p.Nu_2 = @(x) 0.79*(2.06./p.eps_g(x).*p.Re(x).^0.425.*p.Pr(x).^(1/3));
p.Sh = @(x) 1.09./p.eps_g(x).*(p.Re(x).*p.Sc(x)).^(1/3);
p.h_T = @(x) 10;% p.h_T_red*((2.19.*p.Re_bird(x).^(-2/3)).*p.rho_g(x).*p.ug(x).*p.eps_g(x).*p.cp_gas(x)./(p.cp_gas(x)*p.mu/p.k_N2).^(2/3)); % Bird, Re < 50, paragraph 14.5
            %+0.78*p.Re_2(x)^(-.381)
% p.h_T2 = @(x) p.Nu_2(x)*p.k_N2/p.d_p; % Gupta and Thodos, Re > 90
% p.h_T3 = @(x) p.Nu(x)*p.k_N2/p.d_p; % Kramers. Re > 90
p.activ_DR=@(x) min(max(x(:,3)-278,0),1);
p.drying_driving_force = @(x) max(p.Psat_solv(x,p.coeff(x))-x(:,1)./p.MW_solv.*p.MWgas(x).*x(:,5),0);
p.DR = @(x) p.h_M*p.a_V*p.drying_driving_force(x).*p.f(x).*p.activ_DR(x);    

%% Solving routines

p.step_size = p.Z/(p.n_nodes);
x0_1=ones(1,p.n_nodes)*0;
x0_2=ones(1,p.n_nodes)*p.epsL0;
x0_3=ones(1,p.n_nodes)*298;
x0_4=ones(1,p.n_nodes)*298;
p.Pprofile=ones(1,p.n_nodes)*101325;
x0=[x0_1'; x0_2'; x0_3'; x0_4']';

% sparsity matrix - consistent with the FVM upwind differencing scheme
S1=triu(tril(ones(p.n_nodes)),-1); % lower diagonal 
S2=eye(p.n_nodes);
S=repmat(S1,[4,4]);
options = odeset('JPattern',S);

t0=0;
tf=0;
x=[];
t=[];
options_P=optimset('Display','off');

while tf < p.tf   
    tf=tf+step_mom;
    p.Pprofile=fsolve(@momentum_balance,p.Pprofile,options_P,p,x0);
    [t_new,x_new]=ode15s(@model_dryer,t0:5:tf,x0,options,p);
    t0=t_new(end);
    x0=x_new(end,:);
    x=[x;[x_new(1:end-1,:) repmat(p.Pprofile,[size(t_new,1)-1,1])]];
    t=[t;t_new(1:end-1,:)];    
end


%% Graphical output

t_step_graph=graph_interval/5;
z = linspace(p.step_size/2,p.Z-p.step_size/2,p.n_nodes);
legend_label=[];
for i = 1:t_step_graph:length(t)
    legend_label=[legend_label {['t = ' num2str(round(t(i)/60,-1)) ' min']}];
end

xsolv=x(:,1:p.n_nodes);
eps_l=x(:,p.n_nodes+1:2*p.n_nodes);
Tg=x(:,1+2*p.n_nodes:3*p.n_nodes);
Ts=x(:,1+3*p.n_nodes:4*p.n_nodes);
P=x(:,1+4*p.n_nodes:5*p.n_nodes);

figure
plot(z/p.Z,eps_l(1:t_step_graph:end,:)*100,'linewidth',2)
legend(legend_label,'location','best','fontsize',10)
set(gca,'linewidth',1,'fontsize',18,'ytick',0:2.5:15)
axis([0 1 0 12.5])
xlabel([{'Normalized axial position z/Z [-]'}])
ylabel([{'Solvent concentration'} {'in the cake [%]'}])

figure
plot(z,Tg(1:t_step_graph:end,:),'linewidth',1)
legend(legend_label,'location','best','fontsize',10)
set(gca,'linewidth',1,'fontsize',18)
xlabel('Axial position z [m]')
ylabel('Gas phase temperature [K]')

figure
plot(z/p.Z,Ts(1:t_step_graph:end,:)-273.15,'linewidth',2)
legend(legend_label,'location','best','fontsize',10)
set(gca,'linewidth',1,'fontsize',18,'ytick',10:15:85)
axis([0 1 10 85])
xlabel([{'Normalized axial position z/Z [-]'}])
ylabel('Cake temperature [°C]')

figure
plot(z,P(1:t_step_graph:end,:),'linewidth',1)
legend(legend_label,'location','best','fontsize',10)
set(gca,'linewidth',1,'fontsize',18)
xlabel('Axial position z [m]')
ylabel('Pressure [Pa]')

return

% for j = 1 : p.n_nodes
%     for i = 1 :length(t)
% %         xj.(['node_' num2str(j)]).rhog(i)=p.rho_g(x(i,j:p.n_nodes:end));
% %         xj.(['node_' num2str(j)]).states(i,:)=x(i,j:p.n_nodes:end);
% %         xj.(['node_' num2str(j)]).ug(i)=p.ug(x(i,j:p.n_nodes:end));
% %         xj.(['node_' num2str(j)]).Pr(i)=p.Pr(x(i,j:p.n_nodes:end));
% %         xj.(['node_' num2str(j)]).Re(i)=p.Re(x(i,j:p.n_nodes:end));
%         xj.(['node_' num2str(j)]).Re_bird(i)=p.Re_bird(x(i,j:p.n_nodes:end));
% %         xj.(['node_' num2str(j)]).Nu(i)=p.Nu(x(i,j:p.n_nodes:end));
% %         xj.(['node_' num2str(j)]).Nu_2(i)=p.Nu_2(x(i,j:p.n_nodes:end));
%         xj.(['node_' num2str(j)]).h_T1(i)=p.h_T(x(i,j:p.n_nodes:end));
% %         xj.(['node_' num2str(j)]).h_T2(i)=p.h_T2(x(i,j:p.n_nodes:end));
% %         xj.(['node_' num2str(j)]).h_T3(i)=p.h_T3(x(i,j:p.n_nodes:end));
% %         xj.(['node_' num2str(j)]).DR(i)=p.DR(x(i,j:p.n_nodes:end));
%     end
% % % %     figure(1),plot(xj.(['node_' num2str(j)]).rhog),hold on   
% % % %     figure(3),plot(xj.(['node_' num2str(j)]).Pr),hold on
% % % %     figure(5),plot(xj.(['node_' num2str(j)]).Nu),hold on,%plot(xj.(['node_' num2str(j)]).Nu_2),
% % % %     figure(7),plot(xj.(['node_' num2str(j)]).DR),hold on
% end

% j=1;
% figure(1),plot(p.ug(x(:,j:p.n_nodes:end))),ylabel('ug')
% figure(2),plot(xj.(['node_' num2str(j)]).Re_bird),ylabel('Re')
% figure(3),plot(xj.(['node_' num2str(j)]).h_T1),ylabel('h_T')

legend_label=[];
for i = p.step_size/2:p.step_size:p.Z-p.step_size/2;
    legend_label=[legend_label {['z = ' num2str(i) ' m']}];
end
% variables order: |x_solv, eps_L, Tgas, Tsol, P|0...

figure
plot(t,xsolv,'linewidth',1)
% hold on
% plot(z,xsolv_sat(1:t_step_graph:end,:),'linewidth',1)
% legend(legend_label,'location','best','fontsize',10)
set(gca,'linewidth',1,'fontsize',14)
xlabel('Axial position z [m]')
ylabel('Solvent conc. in gas phase [kg/kg]')

figure
plot(t,eps_l,'linewidth',1)
% legend(legend_label,'location','best','fontsize',10)
set(gca,'linewidth',1,'fontsize',14)
xlabel('Axial position z [m]')
ylabel('Liquid phase volum. fraction [-]')

figure
plot(t,Tg,'linewidth',1)
% legend(legend_label,'location','best','fontsize',10)
set(gca,'linewidth',1,'fontsize',14)
xlabel('Axial position z [m]')
ylabel('Gas phase temperature [K]')

figure
plot(t,Ts,'linewidth',1)
% legend(legend_label,'location','best','fontsize',10)
set(gca,'linewidth',1,'fontsize',14)
xlabel('Axial position z [m]')
ylabel('Solid phase temperature [K]')