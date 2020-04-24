function drying_output = model_drying(input,p)
    
    % Grid discretization
    p.number_nodes = p.number_nodes_drying;
    p.nodes_list=1:p.number_nodes;
    p.step_grid_drying=p.L_cake/p.number_nodes;
    p.nodes_drying=linspace(p.step_grid_drying/2,p.L_cake-p.step_grid_drying/2,...
        p.number_nodes);
     
    % Initialization
    p.dPgdz=p.dP_drying/p.L_cake;
%     p.ug = @(x) p.k/p.viscG./(p.eps_g(x))*p.dPgdz;
    p.ug = @(x) p.k/p.viscG*p.dPgdz*sqrt(x(:,3)/298); % added temperature dependence of viscosity
    p.Pprofile=linspace(101325+p.dP_drying-p.dPgdz*p.step_grid_drying/2,...
        101325-p.dPgdz*p.step_grid_drying/2,p.number_nodes);
    S0=input.S_final; 
    
    % Rescale initial saturation on washing grid if deliquoring and washing grids are different
    if abs(input.nodes_deliq(4)-p.nodes_drying(4))>1e-10
        warning('off','all');
        poly_fit_S=polyfit(input.nodes_deliq,S0,4);
        warning('on','all');
        S0=polyval(poly_fit_S,p.nodes_drying);
    end
    
    % Initial conditions
    solv_concG_0=ones(1,p.number_nodes)*0;
    epsL_0=S0*p.E;
    Tg_0=ones(1,p.number_nodes)*298;
%     cp_g_0=p.cp_gas([0,epsL0(1),Tg_0(1),Tg_0(1),1e5]);
    Ts_0=ones(1,p.number_nodes)*298;
    x0=[solv_concG_0'; epsL_0'; Tg_0'; Ts_0']';
    

    % sparsity matrix - consistent with the FVM upwind differencing scheme
%     S1=triu(tril(ones(p.number_nodes)),-1); % lower diagonal
    S1=triu(tril(ones(p.number_nodes),1),-1); % lower diagonal + upper diagonal (HRFVM)
%     S2=eye(p.number_nodes);
    S=repmat(S1,[4,4]);
    options = odeset('JPattern',S);

    t0=0;
    tf=0;
    x=[];
    t=[];
    options_P=optimset('Display','off');

    [t_drying,x]=ode23s(@drying_pde,0:p.time_step_drying:p.t_rot,x0,options,p);
    
    xv_mother_liquor_in_gas=x(:,1:p.number_nodes);
    xv_liquid_phase=x(:,p.number_nodes+1:2*p.number_nodes);
    Tgas=x(:,1+2*p.number_nodes:3*p.number_nodes);
    Tsolid=x(:,1+3*p.number_nodes:4*p.number_nodes);
    
    drying_output.t_drying=t_drying;
    drying_output.xv_mother_liquor_in_gas=xv_mother_liquor_in_gas;
    drying_output.xv_liquid_phase=xv_liquid_phase;
    drying_output.Tgas=Tgas;
    drying_output.Tsolid=Tsolid;
    drying_output.nodes_drying=p.nodes_drying;

end

function dxdt = drying_pde(t,x,p)
    
    x_j=zeros(p.number_nodes,5);
     
    % States matrix
    x=[x; p.Pprofile'];
    for j = 1:p.number_nodes
        x_j(j,:)=x(j:p.number_nodes:end);
    end

    % Discretization: first order upwind scheme
    ug_j=p.ug(x_j);
    rho_g_j=p.rho_g(x_j);
    cp_gas_j=p.cp_gas(x_j);
    DR_j=p.DR(x_j);
    h_T_j=p.h_T(x_j);
    eps_g_j=p.eps_g(x_j);
    
    % fluxes at inlet
    x_jminus1=[0 x_j(1,2) p.Tinlet_drying x_j(1,4) 101325+p.dP_drying]; 
%     x_jminus1=[0 0 p.Tinlet_drying p.Tinlet_drying 101325+p.dP_drying]; 
    flux_eq1_in=0; % the inlet conc of solvent in gas is 0 usually, unless you
                   % have moisture in the air --> x_jminus1(1)*p.ug(x_jminus1)*p.rho_g(x_jminus1);
    flux_eq3_in=x_jminus1(3)*p.cp_gas(x_jminus1);
    
    % fluxes between the nodes
    flux_eq1=x_j(1:p.number_nodes,1);
    flux_eq3=x_j(1:p.number_nodes,3).*cp_gas_j;
%     flux_eq1=fluxes_HRFVM_vL(x_j(1:p.number_nodes,1),p);
%     flux_eq3=fluxes_HRFVM_vL(x_j(1:p.number_nodes,3).*cp_gas_j,p);
    
    % fluxes differences
    dflux_eq1_dz(1)=(flux_eq1(1)-flux_eq1_in)/(p.step_grid_drying);
    dflux_eq3_dz(1)=(flux_eq3(1)-flux_eq3_in)/(p.step_grid_drying);
    dflux_eq1_dz(2:p.number_nodes-1,1)=(flux_eq1(2:p.number_nodes-1)-flux_eq1(1:p.number_nodes-2))/p.step_grid_drying;
    dflux_eq3_dz(2:p.number_nodes-1,1)=(flux_eq3(2:p.number_nodes-1)-flux_eq3(1:p.number_nodes-2))/p.step_grid_drying;  
    % Multiply by 2 if using HRFVM
    dflux_eq1_dz(p.number_nodes,1)=(flux_eq1(p.number_nodes)-flux_eq1(p.number_nodes-1))/p.step_grid_drying; 
    dflux_eq3_dz(p.number_nodes,1)=(flux_eq3(p.number_nodes)-flux_eq3(p.number_nodes-1))/p.step_grid_drying;                

    % mass and energy balances
    
    MB_solv_G= (-ug_j.*rho_g_j.*dflux_eq1_dz+(1-x_j(:,1)).*DR_j)./(rho_g_j.*eps_g_j); % ODE - solvent concentration in G phase
    MB_L=-DR_j/p.rho_liq; % ODE - MB L phase
    EB_G= (-2*rho_g_j.*eps_g_j.*x_j(:,3).*p.dcp_gasdx1(x_j).*MB_solv_G-ug_j.*...
        rho_g_j.*dflux_eq3_dz-h_T_j*p.a_V.*(x_j(:,3)-x_j(:,4))-x_j(:,4).*DR_j)./...
        (rho_g_j.*eps_g_j.*cp_gas_j+rho_g_j.*eps_g_j.*x_j(:,3).*p.dcp_gasdT(x_j)); % ODE - air inlet temperature
    EB_S=(-DR_j*p.latent_heat+h_T_j*p.a_V.*(x_j(:,3)-x_j(:,4))-p.cp_l*p.rho_liq*...
        x_j(:,4).*MB_L)./(p.cp_s*p.rho_sol*(1-p.E)+p.cp_l*p.rho_liq*x_j(:,2)); % ODE - EB L+S
    
    dxdt = [MB_solv_G;MB_L;EB_G;EB_S];  

end

function fluxes=fluxes_HRFVM_vL(x,p)
    r(2:p.number_nodes-1)=(x(2:p.number_nodes-1)-x(1:p.number_nodes-2)+1e-8)./...
        (x(3:p.number_nodes)-x(2:p.number_nodes-1)+1e-8);
    phi_r=((r+abs(r))./(1+abs(r)))'; % the first element is a zero
%     phi_r=[0;(max(max(0,min(1,2*r)),min(r,2)))'];
%     phi_r=[0;(max(0,min(1,r)))'];
    fluxes(1)=0.5*(x(1)+x(2));
    fluxes(2:p.number_nodes-1)=x(2:p.number_nodes-1)+0.5*phi_r(2:p.number_nodes-1).*(x(3:p.number_nodes)-x(2:p.number_nodes-1));
    fluxes(p.number_nodes)=x(end);
end