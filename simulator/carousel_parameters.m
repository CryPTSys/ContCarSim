function p = carousel_parameters(p)
    %% Setup of the components of the liquid phase of the system
    p.names_components=[{'EtOH'}]; % liquid phase components
    p.number_components=length(p.names_components);
    p.number_volatile_components=1; %
    
    %% Liquid phase physical properties    
    % pure components liquid viscosity coefficients (Yaw's) - eq as in visc_liq_components method of carousel_parameters_class_control class
    p.visc_liq_coeff = [-3.1970 7.4084E+02 4.6291E-03 -7.1715E-06] ;         % coefficients ethanol    
    p.MW_components=46.07*1e-3;               % molecular weight ethanol [kg/mol]
    p.rho_liq_components=842;                 % pure ethanol liquid densities [kg/m^3] 
    p.cp_liq_components = 2570;               % pure ethanol specific heat [J/(kg K)]
    p.latent_heat = 846*1e3;                  % pure ethanol latent heat of vaporization [J/kg]

    % Antoine equation coefficients (Perry's) - the eq is coded inside the drying model
    %   Psat = exp(A1+A2/T+A3*log(T)+A4*T^A5;   
    p.coeff_antoine = [74.475, -7164.3,-7.327, 3.1340e-6, 2];
    p.surf_t = 22.39e-3;             % Surface tension [N/m] - EtOH (approximate)
    
    %% Carousel parameters
    p.station_diameter = 0.0152;       % Filter diameter [m]
%     p.Rm = [3e9 3e9 3e9 3e9]';         % Filter medium resistance [m^-1];
                
    %% Gas phase physical properties 
    % assumed same as pure air
    p.k_air = 0.028;                     % air conductivity [W/(m K)]
    p.MW_air = 28.97e-3;                 % air molecular weight [kg/mol]
    cp_air = [29e-3, 0.2199e-5, 0.5723e-8, -2.871e-12]/p.MW_air*1e3; % nitrogen specific heat [J/(kg K)] 
    cp_EtOH= [61.34e-3, 15.72e-5, 8.749e-8, 19.83e-12]/p.MW_components*1e3;
    p.cp_gas_components=[cp_air; cp_EtOH];
    p.visc_gas_phase = 1.8e-5;          % air viscosity at room temperature and 1 atm [Pa s] (temperature effect considered in drying model)
    p.cp=@(T) cp_air(1)+cp_air(2)*(T+273.15)+cp_air(3)*(T+273.15)^2+cp_air(4)*(T+273.15)^3;
    
    %% Solid phase physical properties 
    p.rho_sol = 1293;                % PCM density [kg/m^3]
    p.lambda = 5;                    % pore size index for deliquoring model [-]
    p.cp_s = 2267;                   % Drying - PCM solid specific heat [J/(kg K)] - regressed   
    p.a_V = 126000;                  % Specific surface area [m2/m3] - measured with master sizer
    p.E=0.35;                        % Nominal cake porosity [-]
    p.alpha=2.7e9;                   % Nominal specific cake resistance [m/kg]
    
    %% Kinetics
    p.wl_crit = 0.05;                % Drying - critical impurity content [kgi/kg_cake]
    p.wl_eq = 0.0005;                % Drying - equilibrium impurity content [kgi/kg_cake] 
    p.zeta=199.78;                   % channeling parameter - regressed

    % Filter surface calculation
    p.A = p.station_diameter.^2.*pi./4;  % filtration area [m^2]
    
end