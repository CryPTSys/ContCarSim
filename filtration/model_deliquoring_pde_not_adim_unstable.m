function d = model_deliquoring_pde(cryst_output,filt_output,p)

    %% Deliquoring model calculations - solution obtained solving the PDEs model
    t_deliq=0:.1:p.t_deliq_final;
    L_cake=filt_output.H_cake(end);           % cake height at the end of filtration [m]
    p.k= 1/(p.alpha*p.rho_sol*(1-p.E));      % cake permeability [?]
    p.N_cap_CSD=@(x) (p.E^3*(x).^2*(p.rho_liq*9.81.*L_cake+p.dP))/((1-p.E)^2*L_cake.*p.surf_t);
    p.S_inf=trapz(cryst_output.x,0.155*(1+0.031*p.N_cap_CSD(cryst_output.x).^(-0.49)).*cryst_output.CSD'/p.m0);    % irreducible cake saturation (minimum saturation that can be achieved by displacement of the interstitial liquid by the applied vacuum    
    pb_CSD=@(x) 4.6*(1-p.E)*p.surf_t./(p.E*x);
    p.Pb=trapz(cryst_output.x,pb_CSD(cryst_output.x).*cryst_output.CSD'/p.m0);
    
    % Grid discretization and initial conditions
    p.number_nodes=round(L_cake/p.grid_deliq);
    p.step_grid_deliq=L_cake/p.number_nodes;
    p.nodes_deliq=linspace(p.step_grid_deliq/2,L_cake-p.step_grid_deliq/2,p.number_nodes);
   % dPgdz=p.dP/L_cake;
    p.Pgin=101325;
    p.Pgout=101325-p.dP;
%     p.Pg=linspace(p.Pgin-dPgdz*p.step_grid_deliq/2,...
%         p.Pgout+dPgdz*p.step_grid_deliq/2,length(p.nodes_deliq));
    p.Pg=linspace(p.Pgin,p.Pgout,p.number_nodes+2);
    p.Pg=p.Pg(2:end-1);
    S0=ones(1,p.number_nodes); 
    %Pg0=Pl0+p.Pb*SR0.^(-1/p.lambda);
    
    [t,S_t]=ode15s(@deliquoring_pde,0:.1:p.t_deliq_final,S0,[],p);
    
    S=trapz(p.nodes_deliq,S_t'/(L_cake-p.step_grid_deliq)); % average saturation
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
    d.nodes_deliq=p.nodes_deliq;
   
%     d.SR=SR;
    
end
% HRFVM
function dSdt=deliquoring_pde(t,x,p)
    
    S=x(1:length(p.nodes_deliq));
    SR=(S-p.S_inf)/(1-p.S_inf);    
%     Plout=min(p.Pgout-p.Pb*SR(end)^(-1/p.lambda),p.Pgout);
    ulin=0;       
    Pl=p.Pg'-p.Pb.*SR.^(-1/p.lambda);
%     Plin=min(p.Pgin-p.Pb*SR(1)^(-1/p.lambda),p.Pgin);
    kl=p.k*SR.^((2+3*p.lambda)/p.lambda);
    
%     r_Pl(2:p.number_=(Pl(2:end-1)-Pl(1:end-2)+1e-2)./(Pl(3:end)-Pl(2:end-1)+1e-2);
%     r_Pl_1=(Pl(1)-Plin+1e-6)/(Pl(2)-Pl(1));
%     phi_r_Pl=[(r_Pl_1+abs(r_Pl_1))./(1+abs(r_Pl_1)); (r_Pl+abs(r_Pl))./(1+abs(r_Pl))];
%     phi_r_Pl=[max(0,min(1,r_Pl_1)); max(0,min(1,r_Pl))];
%     fluxPl=Pl(1:end-1)+0.5*phi_r_Pl.*(Pl(2:end)-Pl(1:end-1));
%     dPldz(2:p.number_nodes-1)=(fluxPl(2:end)-fluxPl(1:end-1))/p.step_grid_deliq;
%     dPldz(1)=(fluxPl(1)-Plin)/p.step_grid_deliq;

    flux_Pl_in=Pl(1);
    fluxes_Pl=fluxes_HRFVM_vL(Pl,p);
%     fluxes_Pl(p.number_nodes)=Plout;
    dPldz(1)=(fluxes_Pl(1)-flux_Pl_in)*2/p.step_grid_deliq;
    dPldz(2:p.number_nodes)=(fluxes_Pl(2:p.number_nodes)-fluxes_Pl(1:p.number_nodes-1))/p.step_grid_deliq;
    dPldz(p.number_nodes)=(fluxes_Pl(p.number_nodes)-fluxes_Pl(p.number_nodes-1))*2/p.step_grid_deliq;
    
    ul=-kl./p.visc.*dPldz';
    
%     r_ul=(ul(2:end-1)-ul(1:end-2)+1e-6)./(ul(3:end)-ul(2:end-1)+1e-6);
%     r_ul_1=(ul(1)-ulin+1e-6)/(ul(2)-ul(1));
% %     phi_r_ul=[(r_ul_1+abs(r_ul_1))/(1+abs(r_ul_1)); (r_ul+abs(r_ul))./(1+abs(r_ul))];
%     phi_r_ul=[max(0,min(1,r_ul_1)); max(0,min(1,r_ul))];
%     fluxul=ul(1:end-1)+0.5*phi_r_ul.*(ul(2:end)-ul(1:end-1)); 
%     duldz(2:p.number_nodes-1)=(fluxul(2:end)-fluxul(1:end-1))/p.step_grid_deliq;
%     duldz(1)=(ul(1)-ulin)/p.step_grid_deliq;
%     duldz(p.number_nodes)=(ul(end)-ul(end-1))/p.step_grid_deliq;
%     duldz(2:p.number_nodes)=(ul(2:p.number_nodes)-ul(1:p.number_nodes-1))/p.step_grid_deliq;
%     duldz(1)=(ul(1)-ulin)/(p.step_grid_deliq/2);
%     duldz(2:p.number_nodes)=(ul(2:p.number_nodes)-ul(1:p.number_nodes-1))/p.step_grid_deliq;

    flux_ul_in=ulin;
    fluxes_ul=fluxes_HRFVM_vL(ul,p);
%     fluxes_Pl(p.number_nodes)=Plout;
    duldz(1)=(fluxes_ul(1)-flux_ul_in)/p.step_grid_deliq;
    duldz(2:p.number_nodes)=(fluxes_ul(2:p.number_nodes)-fluxes_ul(1:p.number_nodes-1))/p.step_grid_deliq;
    duldz(p.number_nodes)=(fluxes_ul(p.number_nodes)-fluxes_ul(p.number_nodes-1))*2/p.step_grid_deliq;

    dSdt=-duldz';
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

function fluxes=fluxes_HRFVM_vL(x,p)
    r(2:p.number_nodes-1)=(x(2:p.number_nodes-1)-x(1:p.number_nodes-2)+1e-8)./...
        (x(3:p.number_nodes)-x(2:p.number_nodes-1)+1e-8);
%     phi_r=((r+abs(r))./(1+abs(r)))'; % the first element is a zero
    phi_r=[0;(max(max(0,min(1,2*r)),min(r,2)))'];
    fluxes(1)=0.5*(x(1)+x(2));
    fluxes(2:p.number_nodes-1)=x(2:p.number_nodes-1)+0.5*phi_r(2:p.number_nodes-1).*(x(3:p.number_nodes)-x(2:p.number_nodes-1));
    fluxes(p.number_nodes)=x(end);
end
