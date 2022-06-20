function deliq_output = layer2_model_deliquoring(x_estim,duration_deliq,u)
    
    % deliquoring model used within RTO and for estimating cake saturation
    % at beginning of drying

    % Inputs:
    % x_estim
    % duration_deliq   =   deliquoring duration
    %
    % Output:
    % deliq_output.S = vector of final liquid phase saturation [1 x number_nodes]
    
    %% Deliquoring model calculations - solution obtained solving the PDEs model       
    % Grid discretization variables
    number_nodes=10;
    step_grid_deliq=x_estim.L_cake/number_nodes;
    S0=ones(1,10);

    % Parameters
    visc_liq_coeff = [-3.1970 7.4084E+02 4.6291E-03 -7.1715E-06];
    k=1/(x_estim.alpha*1293*(1-x_estim.E)); 
    S_inf=0.085;
    Pb=4300; 
    
    
    % Calculations    
    visc_liq=10.^(visc_liq_coeff(1)+visc_liq_coeff(2)/295.25+...
                visc_liq_coeff(3)*295.25+visc_liq_coeff(4)*295.25^2)*1e-3; 
    m_solid_initial=x_estim.V_slurry.*x_estim.c_slurry;    % Total solid mass filled into the filter [kg]
    V_solid_initial=m_solid_initial./1293;        % Total solid volume filled into the filter [m^3]
    V_liq_initial=x_estim.V_slurry-V_solid_initial;  % Total liquid volume filled into the filter [m^3]
    V_liquid_pores_end_of_filtration=V_solid_initial*x_estim.E/(1-x_estim.E); % Volume of liquid in the pores at the end of filtration [m3]
    V_filt_final=V_liq_initial-V_liquid_pores_end_of_filtration;   % volume of filtrate at the end of filtration [m3]
    c=(m_solid_initial)/V_filt_final; % Mass of dry cake deposited per unit volume of filtrate [kg sol/ m3 filtrate liquid]more compact way, but same results as Brigi's
    dP_mesh=u.P_compr./(x_estim.alpha*c*V_filt_final/181.4584e-006+x_estim.Rm)*x_estim.Rm; 
    dPgdz=(u.P_compr-dP_mesh)/x_estim.L_cake;
    Pgin=101325+u.P_compr-dPgdz*step_grid_deliq/2; 
    Pgout=101325-dPgdz*step_grid_deliq/2;   
    Pg=linspace(Pgin,Pgout,number_nodes);

    SR0=(S0-S_inf)/(1-S_inf); 
    scaling=1/(step_grid_deliq^4); % an additional scaling factor for pressures, needed for numerical stability
    theta_fin=k*Pb*scaling*duration_deliq/(visc_liq*x_estim.L_cake^2*x_estim.E*(1-S_inf));

    % Sparsity matrix - consistent with the FVM upwind differencing scheme
    SparsMatrix=triu(tril(ones(number_nodes),1),-1); % lower diagonal + upper diagonal
    options = odeset('JPattern',SparsMatrix);
    step_grid_deliq=step_grid_deliq/x_estim.L_cake; %  grid step becomes dimensionless

    % Solver call
    [~,S]=ode15s(@deliquoring_grad_mex,[0 theta_fin],SR0,options,number_nodes,...
       step_grid_deliq, Pb, 5, scaling, Pgin, Pgout, Pg);

%     step_grid_deliq=step_grid_deliq*x_estim.L_cake; % grid step has dimension again
    S=S_inf+S*(1-S_inf); % absolute saturation time profile

    %% Collect outputs in the object deliq_output
    deliq_output.S=S(end,:);    
end
