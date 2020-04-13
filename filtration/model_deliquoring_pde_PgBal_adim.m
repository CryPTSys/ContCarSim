function d = model_deliquoring_pde_PgBal_adim(cryst_output,filt_output,p)

    %% Deliquoring model calculations - solution obtained solving the PDEs model
    t_deliq=0:.1:p.t_deliq_final;
    p.L_cake=filt_output.H_cake(end);           % cake height at the end of filtration [m]
    p.k= 1/(p.alpha*p.rho_sol*(1-p.E));      % cake permeability [?]
    p.N_cap_CSD=@(x) (p.E^3*(x).^2*(p.rho_liq*9.81.*p.L_cake+p.dP))/((1-p.E)^2*p.L_cake.*p.surf_t);
    p.S_inf=trapz(cryst_output.x,0.155*(1+0.031*p.N_cap_CSD(cryst_output.x).^(-0.49)).*cryst_output.CSD'/p.m0);    % irreducible cake saturation (minimum saturation that can be achieved by displacement of the interstitial liquid by the applied vacuum    
    pb_CSD=@(x) 4.6*(1-p.E)*p.surf_t./(p.E*x);
    p.Pb=trapz(cryst_output.x,pb_CSD(cryst_output.x).*cryst_output.CSD'/p.m0);

    % Grid discretization and initial conditions
    p.number_nodes=round(p.L_cake/p.grid_deliq);
    p.step_grid_deliq=p.L_cake/p.number_nodes;
    p.nodes_deliq=linspace(p.step_grid_deliq/2,p.L_cake-p.step_grid_deliq/2,p.number_nodes);
    p.scaling=1/(p.step_grid_deliq^2);
    
    p.Pgin=101325/(p.Pb*p.scaling);
    p.Pgout=p.Pgin-p.dP/(p.Pb*p.scaling);
    p.ulin=0;
    
    dPg0dz=p.dP/p.L_cake/(p.Pb*p.scaling);
    Pg0=linspace(p.Pgin-dPg0dz*p.step_grid_deliq/2,...
        p.Pgout+dPg0dz*p.step_grid_deliq/2,p.number_nodes);
    SR0=ones(1,p.number_nodes)*(0.9999);    
    
    theta_step=p.k*p.Pb*p.scaling*0.1/(p.visc*p.L_cake^2*p.E*(1-p.S_inf));
    theta_fin=p.k*p.Pb*p.scaling*p.t_deliq_final/(p.visc*p.L_cake^2*p.E*(1-p.S_inf));
    
    [~,x]=ode15s(@deliquoring_pde,0:theta_step:theta_fin,[SR0, Pg0.*(1-SR0)],[],p);
    
    SR_t=x(:,1:p.number_nodes);
    SR_t=p.S_inf+SR_t*(1-p.S_inf);
    S=trapz(p.nodes_deliq,SR_t'/(p.L_cake-p.step_grid_deliq)); % average saturation
    
    Pg=x(:,1+p.number_nodes:end)*p.Pb*p.scaling.*(1-x(:,1:p.number_nodes)).^-1;
    V_liq_pores_deliq=S*filt_output.V_liquid_pores_filtration; % Volume of liquid in the pores during deliq as function of time
    V_deliq=(1-S)*filt_output.V_liquid_pores_filtration;   % Volume of filtrate during deliq as function of time
    V_liq_pores_eq=p.S_inf*filt_output.V_liquid_pores_filtration;  % Volume of liquid in the pores at eq   
    Volume_cake=filt_output.V_liquid_pores_filtration+filt_output.m_dry_cake(end)/p.rho_sol;
    solvent_content_vol_deliq=V_liq_pores_deliq/Volume_cake; % Volumetric fraction of solvent content during deliq
    solvent_content_vol_eq=V_liq_pores_eq/Volume_cake; % Volumetric fraction of solvent content at mech eq

    %% Collect outputs in the object d
    d.t_deliq=t_deliq;
    d.V_deliq=V_deliq;
    d.solvent_content_vol_eq=solvent_content_vol_eq;
    d.solvent_content_vol_deliq=solvent_content_vol_deliq;
    d.S=S;
    d.S_inf=p.S_inf;
    d.S_final=SR_t(end,:)*p.E;
    d.nodes_deliq=p.nodes_deliq;
end
% HRFVM
function dYdt=deliquoring_pde(theta,x,p)
    SR=x(1:length(p.nodes_deliq));
    Pg=x(length(p.nodes_deliq)+1:end)./abs(1-SR);

    p.step_grid_deliq=p.step_grid_deliq/p.L_cake; % Dx/L
%     Plin=p.Pgin-p.Pb.*(SR(1)/2).^(-1/p.lambda)/p.Pb/p.scaling;
    
    Pl=Pg-(p.Pb.*SR.^(-1/p.lambda))/p.Pb/p.scaling;
    
    flux_Pl_in=Pl(1);
    fluxes_Pl=fluxes_HRFVM_vL(Pl,p);
    dPldz(1)=(fluxes_Pl(1)-flux_Pl_in)*2/p.step_grid_deliq;
    dPldz(2:p.number_nodes)=(fluxes_Pl(2:p.number_nodes)-fluxes_Pl(1:p.number_nodes-1))/p.step_grid_deliq;
    dPldz(p.number_nodes)=(fluxes_Pl(p.number_nodes)-fluxes_Pl(p.number_nodes-1))*2/p.step_grid_deliq;

    kl=SR.^((2+3*p.lambda)/p.lambda);
    ul=-kl.*dPldz';
    
    flux_ul_in=p.ulin;
    fluxes_ul=fluxes_HRFVM_vL(ul,p);
    duldz(1)=(fluxes_ul(1)-flux_ul_in)/p.step_grid_deliq;
    duldz(2:p.number_nodes)=(fluxes_ul(2:p.number_nodes)-fluxes_ul(1:p.number_nodes-1))/p.step_grid_deliq;
    duldz(p.number_nodes)=(fluxes_ul(p.number_nodes)-fluxes_ul(p.number_nodes-1))*2/p.step_grid_deliq;
    
    flux_Pg_in=Pg(1);
    fluxes_Pg=fluxes_HRFVM_vL(Pg,p);
    dPgdz(1)=(fluxes_Pg(1)-flux_Pg_in)*2/p.step_grid_deliq;
    dPgdz(2:p.number_nodes)=(fluxes_Pg(2:p.number_nodes)-fluxes_Pg(1:p.number_nodes-1))/p.step_grid_deliq;
    dPgdz(p.number_nodes)=(fluxes_Pg(p.number_nodes)-fluxes_Pg(p.number_nodes-1))*2/p.step_grid_deliq;

    kg=(1-SR).^2.*(1-SR.^((2+p.lambda)/p.lambda));
    ug=-kg.*dPgdz';
       
%     Pgug=Pg.*ug;
%     flux_Pgug_in=Pg(1)*ug(1); %approximated
%     fluxes_Pgug=fluxes_HRFVM_vL(Pgug,p);
%     dPgugdz(1)=(fluxes_Pgug(1)-flux_Pgug_in)*2/p.step_grid_deliq;
%     dPgugdz(2:p.number_nodes)=(fluxes_Pgug(2:p.number_nodes)-fluxes_Pgug(1:p.number_nodes-1))/p.step_grid_deliq;
%     dPgugdz(p.number_nodes)=(fluxes_Pgug(p.number_nodes)-fluxes_Pgug(p.number_nodes-1))*2/p.step_grid_deliq;

    flux_ug_in=ug(1); %approximated
    fluxes_ug=fluxes_HRFVM_vL(ug,p);
    dugdz(1)=(fluxes_ug(1)-flux_ug_in)*2/p.step_grid_deliq;
    dugdz(2:p.number_nodes)=(fluxes_ug(2:p.number_nodes)-fluxes_ug(1:p.number_nodes-1))/p.step_grid_deliq;
    dugdz(p.number_nodes)=(fluxes_ug(p.number_nodes)-fluxes_ug(p.number_nodes-1))*2/p.step_grid_deliq;
  
    dSdtheta=-duldz';
%     dPgdtheta=Pg./(1-SR).*dSdtheta-p.visc/p.viscG./(1-SR).*dPgugdz';
    dPgdtheta=-p.visc/p.viscG.*(ug.*dPgdz'+Pg.*dugdz');
    
    dYdt=[dSdtheta;dPgdtheta];

end

function fluxes=fluxes_HRFVM_vL(x,p)
    r(2:p.number_nodes-1)=(x(2:p.number_nodes-1)-x(1:p.number_nodes-2)+eps)./...
        (x(3:p.number_nodes)-x(2:p.number_nodes-1)+eps);
    phi_r=((r+abs(r))./(1+abs(r)))'; % the first element is a zero
%     phi_r=[(max(max(0,min(1,2*r)),min(r,2)))'];
%     phi_r=(max(0,min(2*r,min(2,(1+2*r)/3))))';
    fluxes(1)=(0.5*x(1)+0.5*x(2));
    fluxes(2:p.number_nodes-1)=x(2:p.number_nodes-1)+0.5*phi_r(2:p.number_nodes-1).*(x(3:p.number_nodes)-x(2:p.number_nodes-1));
    fluxes(p.number_nodes)=x(end);
end

% FVM upwind differencing scheme
% function dYdt=deliquoring_pde(t,x,p)
%     S=x(1:length(p.nodes_deliq));
%     Pg=x(length(p.nodes_deliq)+1:end);
%     SR=(S-p.S_inf)/(1-p.S_inf);
% 
%     Pgin=101325;
%     Plin=Pgin-p.Pb*SR(1)^(-1/p.lambda);
%     Pgout=101325-p.dP;
%     ulin=0;
%         
%     Pl=Pg-p.Pb.*SR.^(-1/p.lambda);
%     kl=p.k*SR.^((2+3*p.lambda)/p.lambda);
%     dPldz(1)=(Pl(1)-Plin)/(p.step_grid_deliq/2);
%     dPldz(2:p.number_nodes)=(Pl(2:p.number_nodes)-Pl(1:p.number_nodes-1))/p.step_grid_deliq;
%     ul=-kl./p.visc.*dPldz';
%     duldz(1)=(ul(1)-ulin)/(p.step_grid_deliq/2);
%     duldz(2:p.number_nodes)=(ul(2:p.number_nodes)-ul(1:p.number_nodes-1))/p.step_grid_deliq;
%     kg=p.k*(1-SR).^2.*(1-SR.^((2+p.lambda)/p.lambda));
%     dPgdz(1)=(Pg(1)-Pgin)/(p.step_grid_deliq/2);
%     dPgdz(2:p.number_nodes)=(Pg(2:p.number_nodes)-Pg(1:p.number_nodes-1))/p.step_grid_deliq;  
%     ug=-kg./p.viscG.*dPgdz';
%     dPgugdz(1)=(Pg(1)*ug(1)-Pgin*ug(1))/(p.step_grid_deliq/2);
%     dPgugdz(2:p.number_nodes)=(Pg(2:p.number_nodes).*ug(2:p.number_nodes)-...
%         Pg(1:p.number_nodes-1).*ug(1:p.number_nodes-1))/p.step_grid_deliq; 
%     dSdt=-1./p.E.*duldz';
%     dPgdt=Pg./(1-S)/p.E.*dSdt-1./(1-S)/p.E.*dPgugdz';
%     
%     dYdt=[dSdt;dPgdt];
% 
% end
% 
%     Pl0in
%     Pl0=linspace(101325-dPl0dz*p.step_grid_deliq/2,...
%         101325-p.dP+dPl0dz*p.step_grid_deliq/2,length(p.nodes_deliq));
%     S0=ones(1,length(p.nodes_deliq));
%     SR0=(S0-p.S_inf)/(1-p.S_inf);
%     Pg0=Pl0+p.Pb*SR0.^(-1/p.lambda);