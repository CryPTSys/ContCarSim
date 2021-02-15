function p = carousel_parameters(p)
    %% Setup of the components of the liquid phase of the system
    % default: 3 components =  mother liquor, wash solvent, non volatile impurity
    p.names_components=[{'IPA'},{'EtOH'},{'impurity'}]; % physical properties and state vectors follow this order
    p.number_components=length(p.names_components);
    p.number_volatile_components=2; % the last "p.number_components-p.number_volatile_components" are non-volatile (impurities)
    
    %% Liquid phase physical properties
    
    % pure components liquid viscosity coefficients (Yaw's) - eq as in visc_liq_components method of carousel_parameters_class_control class
    p.visc_liq_coeff=[-7.4681 1.5464e3 1.2131e-2 -1.158e-5  % coefficients component 1
        -3.1970 7.4084E+02 4.6291E-03 -7.1715E-06           % coefficients component 2
        -3.1970 7.4084E+02 4.6291E-03 -7.1715E-06];         % coefficients component 3
    
    p.MW_components=[60.1,46.07,20]'*1e-3;                    % molecular weights [kg/mol] - components 1:3
    p.rho_liq_components=[786, 842, 800]';                 % pure component liquid densities [kg/m^3] - components 1:3  
    p.Di_liq = [1e-9,1e-9,1e-9]';                  % pure compontntmolecular diffusion coefficients [m2/s] - components 1-3
    p.cp_liq_components = [2667 2570 2667]';                % pure component specific heat [J/(kg K)] - components 1-3
    p.latent_heat = [664*1e3 846*1e3 0]';                  % pure component latent heat of vaporization [J/kg] - components 1-3

    % Antoine equation coefficients (Perry's) - the eq is coded inside the drying model
    % Psat = exp(A1+A2/T+A3*log(T)+A4*T^A5;   
    p.coeff_antoine= [76.964, -7623.8, -7.4924, 5.9436e-18, 6,
        74.475, -7164.3,-7.327, 3.1340e-6, 2,
        0, 0, 0, 0, 0];   
    p.surf_t = 22.39e-3;             % Surface tension [N/m] - IPA and EtOH, avg  
    
    %% Carousel parameters
    p.station_diameter = 0.01;       % Filter diameter [m]
%     p.Rm = 2e10;                     % Filter medium resistance [??m^-1];
             
    %% Gas phase physical properties 
    % assumed same as pure air
    p.k_air = 0.028;                    % air conductivity [W/(m K)]
    p.MW_air = 18*1e-3;                 % air molecular weight [kg/mol]
    p.cp_air = [29e-3;0.2199e-5;0.5723e-8;-2.871e-12]/p.MW_air*1e3; % nitrogen specific heat [J/(kg K)] 
    p.visc_gas_phase = 1.8e-5;          % air viscosity at room temperature and 1 atm [Pa s] (temperature effect considered in drying model)

    %% Solid phase physical properties 
    p.rho_sol = 1260;                % PCM density [kg/m^3]
    p.lambda = 5;                    % pore size index for deliquoring model [-]
    p.cp_s = 2200;                   % Drying - (aspirin) solid specific heat [J/(kg K)]   
    p.a_V = 127000;                  % Specific surface area [m2/m3] - measured with master sizer #
    
    %% Kinetics
    p.h_M = [2.25*1e-5/p.a_V 2.25*1e-5/p.a_V 0]'*5;   % Drying kinetic parameter      #
    p.vl_crit = [0.10 0.10 0 ]';                    % Drying - critical impurity content [m3_i/m3]
    p.vl_eq = [5e-5 5e-5 1]';                       % Drying - equilibrium impurity content [m3_l/m3] 
    p.k_ads=0;                                      % Adsorption kinetic parameter
  
    %% Calculations from set operating conditions and parameters
    p.c_inlet=p.conc_from_mass_fr(p.wash_solvent_mass_fr); % wash solvent composition [kg/m3] - components 1-3
    p.A = p.station_diameter.^2.*pi./4;                    % filtration area [m^2]
end

