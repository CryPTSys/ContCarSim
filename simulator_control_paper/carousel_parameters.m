function p = carousel_parameters(p)
    %% Setup of the components of the liquid phase of the system
    % default: 3 components =  mother liquor, wash solvent, non volatile impurity
    p.names_components=[{'EtOH'}]; % physical properties and state vectors follow this order
    p.number_components=length(p.names_components);
    p.number_volatile_components=1; % the last "p.number_components-p.number_volatile_components" are non-volatile (impurities)
    
    %% Liquid phase physical properties
    
    % pure components liquid viscosity coefficients (Yaw's) - eq as in visc_liq_components method of carousel_parameters_class_control class
    p.visc_liq_coeff=[-3.1970 7.4084E+02 4.6291E-03 -7.1715E-06 ] ;         % coefficients component 1
    
    p.MW_components=46.07*1e-3;                    % molecular weights [kg/mol] - components 1:3
    p.rho_liq_components=842;                 % pure component liquid densities [kg/m^3] - components 1:3  
    p.Di_liq = 1e-9;                  % pure compontntmolecular diffusion coefficients [m2/s] - components 1-3
    p.cp_liq_components = 2570;                % pure component specific heat [J/(kg K)] - components 1-3
    p.latent_heat = 846*1e3;                  % pure component latent heat of vaporization [J/kg] - components 1-3

    % Antoine equation coefficients (Perry's) - the eq is coded inside the drying model
    % Psat = exp(A1+A2/T+A3*log(T)+A4*T^A5;   
    p.coeff_antoine= [74.475, -7164.3,-7.327, 3.1340e-6, 2];
    p.surf_t = 22.39e-3;             % Surface tension [N/m] - IPA and EtOH, avg  
    
    %% Carousel parameters
    p.station_diameter = 0.0152;       % Filter diameter [m]
    p.Rm = [3e9 3e9 3e9 3e9]';                     % Filter medium resistance [??m^-1];
                
    %% Gas phase physical properties 
    % assumed same as pure air
    p.k_air = 0.028;                    % air conductivity [W/(m K)]
    p.MW_air = 28.97e-3;                 % air molecular weight [kg/mol]
    cp_air = [29e-3, 0.2199e-5, 0.5723e-8, -2.871e-12]/p.MW_air*1e3; % nitrogen specific heat [J/(kg K)] 
    cp_EtOH= [61.34e-3, 15.72e-5, 8.749e-8, 19.83e-12]/p.MW_components*1e3;
    p.cp_gas_components=[cp_air; cp_EtOH];
    p.visc_gas_phase = 1.8e-5;          % air viscosity at room temperature and 1 atm [Pa s] (temperature effect considered in drying model)
    p.cp=@(T) cp_air(1)+cp_air(2)*(T+273.15)+cp_air(3)*(T+273.15)^2+cp_air(4)*(T+273.15)^3;
    
    %% Solid phase physical properties 
    p.rho_sol = 1293;                % PCM density [kg/m^3]
    p.lambda = 5;                    % pore size index for deliquoring model [-]
    p.cp_s = 2267;                   % Drying - (aspirin) solid specific heat [J/(kg K)]   
    p.a_V = 126000;                  % Specific surface area [m2/m3] - measured with master sizer #
    
    %% Kinetics
%     p.h_M = 2e-9;   % Drying kinetic parameter      #
    p.vl_crit = 0.05;                    % Drying - critical impurity content [m3_i/m3]
    p.vl_eq = 0.0005;                       % Drying - equilibrium impurity content [m3_l/m3] 
    p.k_ads=0;                                      % Adsorption kinetic parameter
%     p.h_T=0.7;
    p.zeta=199.78;
  
    %% Calculations from set operating conditions and parameters
%     p.c_inlet=p.conc_from_mass_fr(p.wash_solvent_mass_fr); % wash solvent composition [kg/m3] - components 1-3
    p.A = p.station_diameter.^2.*pi./4;                    % filtration area [m^2]
end