function dxdt = model_dryer_FVM_noP_vectorized(t,x,p)
    
    x_j=zeros(p.n_nodes,5);
     
    % States matrix
    x=[x; p.Pprofile'];
    for j = 1:p.n_nodes
        x_j(j,:)=x(j:p.n_nodes:end);
    end

    % Discretization: first order upwind scheme
    ug_j=p.ug(x_j);
    rho_g_j=p.rho_g(x_j);
    cp_gas_j=p.cp_gas(x_j);
    DR_j=p.DR(x_j);
    h_T_j=p.h_T(x_j);
    eps_g_j=p.eps_g(x_j);
    
    % fluxes at inlet
    x_jminus1=[0 x_j(1,2) p.Tg_inlet x_j(1,4) p.Pprofile(1)]; 
    flux_eq1_in=x_jminus1(1)*p.ug(x_jminus1)*p.rho_g(x_jminus1);
    flux_eq3_in=x_jminus1(3)*p.ug(x_jminus1)*p.rho_g(x_jminus1)*p.cp_gas(x_jminus1);
    
    % fluxes between the nodes
    flux_eq1=x_j(1:p.n_nodes,1).*ug_j.*rho_g_j;
    flux_eq3=x_j(1:p.n_nodes,3).*ug_j.*rho_g_j.*cp_gas_j;
    
    % fluxes differences
    diff_flux1(1)=(flux_eq1(1)-flux_eq1_in)/(p.step_size/2);
    diff_flux3(1)=(flux_eq3(1)-flux_eq3_in)/(p.step_size/2);
    diff_flux1(2:p.n_nodes,1)=(flux_eq1(2:p.n_nodes)-flux_eq1(1:p.n_nodes-1))/p.step_size;
    diff_flux3(2:p.n_nodes,1)=(flux_eq3(2:p.n_nodes)-flux_eq3(1:p.n_nodes-1))/p.step_size;                
  
    % mass and energy balances
    
    MB_solv_G= (-diff_flux1+DR_j)./(rho_g_j.*eps_g_j); % ODE - solvent concentration in G phase
    MB_L=-DR_j/p.rho_l; % ODE - MB L phase
    EB_G= (- diff_flux3-h_T_j*p.a_V.*(x_j(:,3)-x_j(:,4)))./(rho_g_j.*eps_g_j.*cp_gas_j); % ODE - air inlet temperature
    EB_S=(-DR_j*p.lambda+h_T_j*p.a_V.*(x_j(:,3)-x_j(:,4)))./(p.cp_s*...
        p.rho_s*p.eps_s+p.cp_l*p.rho_l*x_j(:,2)); % ODE - EB L+S
   
    dxdt = [MB_solv_G;MB_L;EB_G;EB_S];  

end