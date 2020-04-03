%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulator for the convective dryer model %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc

% improve: hM, hT
%% Parameters
% simulation inputs
p.n_nodes = 35;
p.tf = 1000; % process duration [s]
graph_interval=100;
step_mom = 100; % time step for recalculating the pressure profile
p.Z = 0.05; % dryer length [m]
p.step_size = p.Z/(p.n_nodes);
p.D_dryer = 0.01; % dryer diameter [m]
p.epsL0 = 0.10; % volumetric fraction

% fixed parameters
p.eps_l_eq = 0; % equilibrium solvent content [m3_l/m3]
p.rho_s = 1400; % solid density [kg_s/m3] - aspirin
% p.rho_s = 1300; % solid density [kg_s/m3] - paracetamol
p.A_dryer = pi/4*p.D_dryer^2; % dryer cross-section [m2]
% p.cp_s = 190/150*1e3; % solid specific heat [J/(kg K)] paracetamol
p.cp_s = 2200; % aspirin % solid specific heat [J/(kg K)]
%p.h_T = 30; % heat transfer coefficient [1/s?] - model it!
% A = 8.07131; % Antoine constant #1
% B = 1730.63; % Antoine constant #2
% C = 233.426; % Antoine constant #3
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
p.a_V = 3500; %6*(1-p.eps_s)/p.d_p/5; % specific area [m2/m3];


x0_1=ones(1,p.n_nodes)*0;
x0_2=ones(1,p.n_nodes)*p.epsL0;
x0_3=ones(1,p.n_nodes)*298;
x0_4=ones(1,p.n_nodes)*298;
p.Pprofile=ones(1,p.n_nodes)*101325;
x0=[x0_1'; x0_2'; x0_3'; x0_4']';

% mass matrix
% M=repmat([1 1 1 1],[1 p.n_nodes]); % with new inlet BCs (on the inlet)
% % M([1,2*p.n_nodes+1])=0; % with old inlet BCs (on the first node)
% M=diag(M);
% options = odeset('Mass',M); % not used, they are all ODEs\

% sparsity matrix - consistent with the FVM upwind differencing scheme
S1=triu(tril(ones(p.n_nodes)),-1); % lower diagonal 
S2=eye(p.n_nodes);
% S=repmat([S1 S2 S1 S2],[4,1]);
S=repmat(S1,[4,4]);
% S=sparse(S); % convert in sparse form the sparsity pattern
options = odeset('JPattern',S);
options_P=optimset('Display','off');

% fixed methods
p.switch_antoine=@(x) 1-min(max(x(:,3)-60-273.15,0),1);
p.coeff= @(x) [p.A1,p.B1,p.C1].*p.switch_antoine(x)+[p.A2,p.B2,p.C2].*(1-p.switch_antoine(x)); 
p.cp = @(x,coeff) (coeff(1)+coeff(2)*x(:,3)+coeff(3)*x(:,3).^2+coeff(4)*x(:,3).^3); % [ J/kg K]
p.MWgas = @(x) (1-x(:,1))*p.MW_N2+x(:,1)*p.MW_solv;
p.rho_g=@(x) x(:,5)./(8.314*x(:,3)).*p.MWgas(x); % gas phase density [kg/m3]
p.cp_gas = @(x) p.cp(x,p.cp_N2).*(1-x(:,1))+p.cp(x,p.cp_solventG).*x(:,1);
p.Psat_solv = @(x,coeff) 10.^(coeff(:,1)-coeff(:,2)./(coeff(:,3)+(x(:,3)-273.15)))*133.322; % saturation pressure [Pa] 
p.activ_DR=@(x) min(max(x(:,3)-278,0),1);

%% Parameter sensitivity routines

% definition of intervals of variability of the parameters

param_sens;
sens.length=[];

for i = 1:length(fields(sens))-2
    sens.length.(cell2str(sens.names(i)))=length(sens.(cell2str(sens.names(i))));
end

% calculation of drying time
for i = 1:length(fields(sens))-2
    for  j = 1 :sens.length.(cell2str(sens.names(i)))
        % adjustable parameters: reference values in script p_ref
        p_ref;
           
        p.(cell2str(sens.names(i)))=sens.(cell2str(sens.names(i)))(j);
       
        % methods depending on the values of the adjustable parameters
        p.eps_g = @(x) 1-p.eps_s-x(:,2);
        p.f=@(x) min((x(:,2)-p.eps_l_eq)/(p.eps_l_crit-p.eps_l_eq),1); % drying rate efficiency function [adim]
        p.Vg = @(x) p.Vg_N * (101325./x(:,5)).*(x(:,3)/273.15); % flowrate of drying nitrogen in [m3/s]
        p.ug = @(x) p.Vg(x)/p.A_dryer./p.eps_g(x); % gas velocity [m/s]
        p.drying_driving_force = @(x) max(p.Psat_solv(x,p.coeff(x))-x(:,1)./p.MW_solv.*p.MWgas(x).*x(:,5),0);
        p.DR = @(x) p.h_M*p.a_V*p.drying_driving_force(x).*p.f(x).*p.activ_DR(x);    
        p.h_T = @(x) p.h_T_value;
        % Solving routines
        t0=0;
        x0=[x0_1'; x0_2'; x0_3'; x0_4']';
        tf=0;
        x=[];
        t=[];
        flag = 0;
        while flag == 0  
            tf=tf+step_mom;
            p.Pprofile=fsolve(@MomBal_vectorized,p.Pprofile,options_P,p,x0);
        %     p.Pprofile=repmat(1e6,[1,p.n_nodes]);
            [t_new,x_new]=ode15s(@model_dryer_FVM_noP_vectorized,t0:5:tf,x0,options,p);
            t0=t_new(end);
            x0=x_new(end,:);
            x=[x;[x_new(1:end-1,:) repmat(p.Pprofile,[size(t_new,1)-1,1])]];
            eps_l_new=x_new(:,p.n_nodes+1:2*p.n_nodes);
            time_steps=1:length(t_new);
            max_moist_per_time_inst=max(eps_l_new');
            time_steps=time_steps(max_moist_per_time_inst<0.005);
            flag = min(0,[(max_moist_per_time_inst(min(time_steps))-0.005)' length(t_new+1)]);
        end
        results.(cell2str(sens.names(i))).(['value_' num2str(j)]).t=0:5:p.tf;
        results.(cell2str(sens.names(i))).(['value_' num2str(j)]).x = x;
        results.(cell2str(sens.names(i))).(['value_' num2str(j)]).drying_time=t_new(min(time_steps)); % seconds
    end 
end
save(['results_3'],'results')
 %% Graphical output
% 
% for j = 1 : p.n_nodes
%     for i = 1 :length(t)
% %         xj.(['node_' num2str(j)]).rhog(i)=p.rho_g(x(i,j:p.n_nodes:end));
% %         xj.(['node_' num2str(j)]).states(i,:)=x(i,j:p.n_nodes:end);
%         xj.(['node_' num2str(j)]).ug(i)=p.ug(x(i,j:p.n_nodes:end));
% %         xj.(['node_' num2str(j)]).Pr(i)=p.Pr(x(i,j:p.n_nodes:end));
% %         xj.(['node_' num2str(j)]).Re(i)=p.Re(x(i,j:p.n_nodes:end));
% %         xj.(['node_' num2str(j)]).Re_bird(i)=p.Re_bird(x(i,j:p.n_nodes:end));
% %         xj.(['node_' num2str(j)]).Nu(i)=p.Nu(x(i,j:p.n_nodes:end));
% %         xj.(['node_' num2str(j)]).Nu_2(i)=p.Nu_2(x(i,j:p.n_nodes:end));
% %         xj.(['node_' num2str(j)]).h_T1(i)=p.h_T(x(i,j:p.n_nodes:end));
% %         xj.(['node_' num2str(j)]).h_T2(i)=p.h_T2(x(i,j:p.n_nodes:end));
% %         xj.(['node_' num2str(j)]).h_T3(i)=p.h_T3(x(i,j:p.n_nodes:end));
% %         xj.(['node_' num2str(j)]).DR(i)=p.DR(x(i,j:p.n_nodes:end));
%     end
% % %     figure(1),plot(xj.(['node_' num2str(j)]).rhog),hold on   
% % %     figure(3),plot(xj.(['node_' num2str(j)]).Pr),hold on
% % %     figure(5),plot(xj.(['node_' num2str(j)]).Nu),hold on,%plot(xj.(['node_' num2str(j)]).Nu_2),
% % %     figure(7),plot(xj.(['node_' num2str(j)]).DR),hold on
% end
% 
% j=1;
% figure(1),plot(p.ug(x(:,j:p.n_nodes:end))),ylabel('ug')
% % figure(2),plot(xj.(['node_' num2str(j)]).Re_bird),ylabel('Re')
% % figure(3),plot(xj.(['node_' num2str(j)]).h_T1),ylabel('h_T')
% return
% t_step_graph=graph_interval/5;
% z = linspace(p.step_size/2,p.Z-p.step_size/2,p.n_nodes);
% legend_label=[];
% for i = 1:t_step_graph:length(t)
%     legend_label=[legend_label {['t = ' num2str(t(i)) ' s']}];
% end
% 
% xsolv=x(:,1:p.n_nodes);
% eps_l=x(:,p.n_nodes+1:2*p.n_nodes);
% Tg=x(:,1+2*p.n_nodes:3*p.n_nodes);
% Ts=x(:,1+3*p.n_nodes:4*p.n_nodes);
% P=x(:,1+4*p.n_nodes:5*p.n_nodes);
% 
% % for j = 1:p.n_nodes
% %     for tt = 1:length(t)
% %         x_t=[xsolv(tt,j),eps_l(tt,j),Tg(tt,j),Ts(tt,j),P(tt,j)];
% %         xsolv_sat(tt,j)=p.MW_solv/p.MWgas(x_t)/P(tt,j)*p.Psat_solv(x_t); 
% %     end
% % end
% 
% figure
% plot(z,xsolv(1:t_step_graph:end,:),'linewidth',1)
% % hold on
% % plot(z,xsolv_sat(1:t_step_graph:end,:),'linewidth',1)
% legend(legend_label,'location','best','fontsize',10)
% set(gca,'linewidth',1,'fontsize',14)
% xlabel('Axial position z [m]')
% ylabel('Solvent conc. in gas phase [kg/kg]')
% 
% figure
% plot(z,eps_l(1:t_step_graph:end,:),'linewidth',1)
% legend(legend_label,'location','best','fontsize',10)
% set(gca,'linewidth',1,'fontsize',14)
% xlabel('Axial position z [m]')
% ylabel('Liquid phase volum. fraction [-]')
% 
% figure
% plot(z,Tg(1:t_step_graph:end,:),'linewidth',1)
% legend(legend_label,'location','best','fontsize',10)
% set(gca,'linewidth',1,'fontsize',14)
% xlabel('Axial position z [m]')
% ylabel('Gas phase temperature [K]')
% 
% figure
% plot(z,Ts(1:t_step_graph:end,:),'linewidth',1)
% legend(legend_label,'location','best','fontsize',10)
% set(gca,'linewidth',1,'fontsize',14)
% xlabel('Axial position z [m]')
% ylabel('Solid phase temperature [K]')
% 
% figure
% plot(z,P(1:t_step_graph:end,:),'linewidth',1)
% legend(legend_label,'location','best','fontsize',10)
% set(gca,'linewidth',1,'fontsize',14)
% xlabel('Axial position z [m]')
% ylabel('Pressure [Pa]')
% 
% return
% 
% legend_label=[];
% for i = p.step_size/2:p.step_size:p.Z-p.step_size/2;
%     legend_label=[legend_label {['z = ' num2str(i) ' m']}];
% end
% % variables order: |x_solv, eps_L, Tgas, Tsol, P|0...
% 
% figure
% plot(t,xsolv,'linewidth',1)
% % hold on
% % plot(z,xsolv_sat(1:t_step_graph:end,:),'linewidth',1)
% % legend(legend_label,'location','best','fontsize',10)
% set(gca,'linewidth',1,'fontsize',14)
% xlabel('Axial position z [m]')
% ylabel('Solvent conc. in gas phase [kg/kg]')
% 
% figure
% plot(t,eps_l,'linewidth',1)
% % legend(legend_label,'location','best','fontsize',10)
% set(gca,'linewidth',1,'fontsize',14)
% xlabel('Axial position z [m]')
% ylabel('Liquid phase volum. fraction [-]')
% 
% figure
% plot(t,Tg,'linewidth',1)
% % legend(legend_label,'location','best','fontsize',10)
% set(gca,'linewidth',1,'fontsize',14)
% xlabel('Axial position z [m]')
% ylabel('Gas phase temperature [K]')
% 
% figure
% plot(t,Ts,'linewidth',1)
% % legend(legend_label,'location','best','fontsize',10)
% set(gca,'linewidth',1,'fontsize',14)
% xlabel('Axial position z [m]')
% ylabel('Solid phase temperature [K]')

%% Results extraction
for w = 1:length(fields(sens))-2
    dr_vector=[];
    for i = 1:sens.length.(cell2str(sens.names(w)))
        names=[cell2str(sens.names(w)) {['value_' num2str(i)]} {'drying_time'}];
        dr_vector(i)=results.(cell2str(names(1))).(cell2str(names(2))).(cell2str(names(3)));
    end
    figure
    plot(sens.(cell2str(names(1))),dr_vector/60,'linewidth',1.5),xlabel(cell2str(sens.names(w))),ylabel('Drying time [min]')
    set(gca,'fontsize',20,'linewidth',1)
end