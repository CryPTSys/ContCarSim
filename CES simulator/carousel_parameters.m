function p = carousel_parameters(cryst_output,p)
    %% System components setup

    names_components=[{'IPA'},{'EtOH'},{'impurity'}]; % mother liquor, wash solvent, impurities in mother liquor
    p.number_components=length(names_components);
    p.number_volatile_components=2;
    p.visc_liq_coeff=[-7.4681 1.5464e3 1.2131e-2 -1.158e-5 % from Yaws
        -3.1970 7.4084E+02 4.6291E-03 -7.1715E-06
        -3.1970 7.4084E+02 4.6291E-03 -7.1715E-06];
    p.MW_components=[60.1,46.07,20]'*1e-3; % kg/mol
    p.rho_liq_components=[786, 842, 800]'; % Liquid densities |mother liquor,washing solvent| [kg/m^3]
    % Washing
    p.c_inlet=[0 p.rho_liq_components(2) 0]';  % Washing solvent composition [kg_i/m3_liq]
    p.Di_liq = [1e-9,1e-9,1e-9]';             % Washing - diffusivities in liquid phase of volatile species [m2/s]
    % Drying
%     p.Di_gas = [0.4e-4 0.4e-4 0]'; % diffusivity of species into air [m2/s]
    p.cp_liq_components = [2667 2570 2667]';   % liquid specific heat [J/(kg K)]
    % we assume cp_gas = cp_air; p.cp_gas_components = [33.46e-3;0.688e-5;0.7604e-8;-3.593e-12]./p.MW_components*1e3; % solvent (steam) specific heat in the gas phase [J/(kg K)] 
    p.vl_crit = [0.1 0.1 0 ]';             %@ Drying - critical impurity content [m3_i/m3]
    p.vl_eq = [1e-12 1e-12 1]';                  %@ Drying - equilibrium impurity content [m3_l/m3] 
    p.latent_heat = [664*1e3 846*1e3 0]'; % solvent latent heat of vaporization J/kg] - model it
    
    % Antoine equation coefficients (Perry's) - the eq is coded inside the drying model
    % Psat = exp(A1+A2/T+A3*log(T)+A4*T^A5;   
    p.coeff_antoine= [76.964, -7623.8, -7.4924, 5.9436e-18, 6,
        74.475, -7164.3,-7.327, 3.1340e-6, 2,
        0, 0, 0, 0, 0];
    

 %% Filtration/deliquoring/washing parameters
    p.station_diameter = 0.01;       % Filter diameter [m]
    p.Rm = 3e10;                   % Filter medium resistance [??m^-1];
    p.rho_sol = 1260;                % Crystal density [kg/m^3]
    p.visc_gas_phase = 1.8e-5;       % Air viscosity at room temperature and 1 atm [Pa s]
    p.surf_t = 22.39e-3;             % Surface tension [N/m]  - ethanol
    p.lambda = 5;                    % Deliquoring - pore size index for deliquoring [-]
    p.k0 = 0.15;                     % Washing - diffusive backflux parameter
    p.gamma = 5;                     % Washing - diffusive back-flux parameter
    
    %% Calculation of cake properties 
    p.A = p.station_diameter.^2.*pi./4;                        % Filtration area [m^2]
    p.E=porosity_function_shape(p.station_diameter,cryst_output.CSD,cryst_output.x,1);           
    p.m0=trapz(cryst_output.x,cryst_output.CSD); 
    p.m1=trapz(cryst_output.x,cryst_output.CSD.*cryst_output.x); 
    p.m2=trapz(cryst_output.x,cryst_output.CSD.*(cryst_output.x.^2)); 
    p.m3=trapz(cryst_output.x,cryst_output.CSD.*(cryst_output.x.^3)); 
    p.alpha=sum(p.alpha_CSD(cryst_output.x).*cryst_output.CSD);

    p.k=1/(p.alpha*p.rho_sol*(1-p.E));      % cake permeability [?]
    p.a_V=6*p.m2/p.m3; % 109000; %
    p.x=cryst_output.x;
    p.CSD=cryst_output.CSD;
    
    %% Drying
    % Drying parameters
    % Gas phase parameters
    p.k_air = 0.028; % air conductivity [W/(m K)] IT'S AIR EVEN IF NAME IS air!
    p.MW_air = 18*1e-3; % air molecular weight [kg/mol] IT'S AIR EVEN IF NAME IS air!
    p.cp_air = [29e-3;0.2199e-5;0.5723e-8;-2.871e-12]/p.MW_air*1e3; % nitrogen specific heat [J/(kg K)] 
%     p.h_T=1000; % not needed anymore since when we use only one EB
    p.h_M = [2.25*1e-5/p.a_V 2.25*1e-5/p.a_V 0]'*1; %             %@ Drying - mass transfer coefficient
    
    % Solid phase properties
    p.cp_s = 1191;                   % Drying - (paracetamol) solid specific heat [J/(kg K)]
    p.Di_gas = {[-0.08697, 4.7538E-04, 5.0525E-07], [-0.10107, 5.6275E-4, 5.8314E-7],...
        [-0.10107, 5.6275E-4, 5.8314E-7]}; % diffusivity in air of volatile species (Yaws)

end

