function d = model_deliquoring_pde_adim(cryst_output,filt_output,p)

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
    p.nodes_list=1:p.number_nodes;
    p.step_grid_deliq=p.L_cake/p.number_nodes;
    p.nodes_deliq=linspace(p.step_grid_deliq/2,p.L_cake-p.step_grid_deliq/2,p.number_nodes);
   
    dPgdz=p.dP/p.L_cake;
    p.Pgin=101325;
    p.Pgout=101325-p.dP;   
    p.Pg=linspace(p.Pgin-dPgdz*p.step_grid_deliq/2,...
        p.Pgout+dPgdz*p.step_grid_deliq/2,p.number_nodes);
    SR0=ones(1,p.number_nodes); 
    p.scaling=1/(p.step_grid_deliq^2); % an additional scaling factor for pressures
    theta_step=p.k*p.Pb*p.scaling*0.1/(p.visc*p.L_cake^2*p.E*(1-p.S_inf));
    theta_fin=p.k*p.Pb*p.scaling*p.t_deliq_final/(p.visc*p.L_cake^2*p.E*(1-p.S_inf));
    
    % sparsity matrix - consistent with the FVM upwind differencing scheme
%     SparsMatrix=tril(ones(p.number_nodes),1); % lower diagonal 
    options = [];%odeset('JPattern',SparsMatrix);
    [~,S_t]=ode15s(@deliquoring_pde,0:theta_step:theta_fin,SR0,options,p);
    
    S_t=p.S_inf+S_t*(1-p.S_inf);
    S=trapz(p.nodes_deliq,S_t'/(p.L_cake-p.step_grid_deliq)); % average saturation
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
    d.S_final=S_t(end,:)*p.E;
    d.nodes_deliq=p.nodes_deliq;
end

% HRFVM
function dSdtheta=deliquoring_pde(~,x,p)
    p.step_grid_deliq=p.step_grid_deliq/p.L_cake;
    
    SR=x(1:length(p.nodes_deliq));
    Plin=(p.Pgin-p.Pb.*SR(1).^(-1/p.lambda))/p.Pb/p.scaling;
    Plout=(min(p.Pgout-p.Pb*SR(end)^(-1/p.lambda),p.Pgout))/p.Pb/p.scaling;
    ulin=0; 
    
    % Neglect transient
%     Pl=(p.Pg'-p.Pb.*SR.^(-1/p.lambda))/p.Pb/p.scaling;
    
    % Approximate transient before air breakthrough
    if SR(end)<1-1e-6
        Pl=(p.Pg'-p.Pb.*SR.^(-1/p.lambda))/p.Pb/p.scaling;
    else
        Pl=linspace(p.Pg(1)-p.Pb.*SR(1).^(-1/p.lambda),p.Pg(end),p.number_nodes)/p.Pb/p.scaling;
        Pl=Pl';
    end

    Pl=Pl(:);
    kl=SR.^((2+3*p.lambda)/p.lambda);
    
    % HRFVM
    flux_Pl_in=Pl(1);
    fluxes_Pl=fluxes_HRFVM_vL(Pl,p);
    dPldz(1)=(fluxes_Pl(1)-flux_Pl_in)*2/p.step_grid_deliq;
    dPldz(2:p.number_nodes)=(fluxes_Pl(2:p.number_nodes)-fluxes_Pl(1:p.number_nodes-1))/p.step_grid_deliq;
    dPldz(p.number_nodes)=(fluxes_Pl(p.number_nodes)-fluxes_Pl(p.number_nodes-1))*2/p.step_grid_deliq;
    
    % FVM - doesn't work, to be refined
%     flux_Pl_in=Plin;
%     fluxes_Pl=fluxes_FVM_upwind(Pl,p);
%     fluxes_Pl(1)=0.5*(Pl(1)+Plin);
%     dPldz(1)=(fluxes_Pl(1)-flux_Pl_in)*2/p.step_grid_deliq;
%     dPldz(2:p.number_nodes)=(fluxes_Pl(2:p.number_nodes)-fluxes_Pl(1:p.number_nodes-1))/p.step_grid_deliq;
%     dPldz(p.number_nodes)=(fluxes_Pl(p.number_nodes)-fluxes_Pl(p.number_nodes-1))*2/p.step_grid_deliq;

    ul=-kl.*dPldz';
    
    % HRFVM
    flux_ul_in=ulin;
    fluxes_ul=fluxes_HRFVM_vL(ul,p);
    duldz(1)=(fluxes_ul(1)-flux_ul_in)/p.step_grid_deliq;
    duldz(2:p.number_nodes)=(fluxes_ul(2:p.number_nodes)-fluxes_ul(1:p.number_nodes-1))/p.step_grid_deliq;
    duldz(p.number_nodes)=(fluxes_ul(p.number_nodes)-fluxes_ul(p.number_nodes-1))*2/p.step_grid_deliq;

    % FVM - doesn't work
%     flux_ul_in=ulin;
%     fluxes_ul=fluxes_FVM_upwind(ul,p);
%     fluxes_ul(1)=ul(1);
%     duldz(1)=(fluxes_ul(1)-flux_ul_in)*2/p.step_grid_deliq;
%     duldz(2:p.number_nodes)=(fluxes_ul(2:p.number_nodes)-fluxes_ul(1:p.number_nodes-1))/p.step_grid_deliq;
% %     duldz(p.number_nodes)=(fluxes_ul(p.number_nodes)-fluxes_ul(p.number_nodes-1))*2/p.step_grid_deliq;

    dSdtheta=-duldz';
end

function fluxes=fluxes_HRFVM_vL(x,p)
    r(2:p.number_nodes-1)=(x(2:p.number_nodes-1)-x(1:p.number_nodes-2)+1e-8)./...
        (x(3:p.number_nodes)-x(2:p.number_nodes-1)+1e-8);
    phi_r=((r+abs(r))./(1+abs(r)))'; % the first element is a zero
%     phi_r=[0;(max(max(0,min(1,2*r)),min(r,2)))'];
    fluxes(1)=0.5*(x(1)+x(2));
    fluxes(2:p.number_nodes-1)=x(2:p.number_nodes-1)+0.5*phi_r(2:p.number_nodes-1).*(x(3:p.number_nodes)-x(2:p.number_nodes-1));
    fluxes(p.number_nodes)=x(end);
end

function fluxes=fluxes_FVM_upwind(x,p)
    fluxes(1)=0.5*(x(1)+x(2));
    fluxes(2:p.number_nodes)=x(1:p.number_nodes-1);
end
