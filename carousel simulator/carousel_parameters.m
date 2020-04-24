function p = carousel_parameters(cryst_output,p)
    %% System components setup
    % in this version we neglect the effect on the gas phase cp of the vaporized species
    p.names_components=[{'H20'},{'EtOH'},{'impurity'}]; % mother liquor, washing solvent, impurities in mother liquor
    p.visc_liq_components = [0.89e-3 1E-03 1E-03];  % Liquid viscosity [Pa s]
    p.MW_components=[18,46.07,20]*1e-3; % kg/mol
    p.rho_liq_components=[1000, 842, 900]; % Liquid densities |mother liquor,washing solvent| [kg/m^3]
    % Washing
    p.c_inlet=[0 p.rho_liq_components(2) 0];  % Washing solvent composition [kg_i/m3_liq]
    p.Di_liq = [1.28e-9,1.28e-9,1.28e-9];             % Washing - diffusivities volatile species [m2/s]
    % Drying
    p.Di_gas = [0.4e-4 0.4e-4 0]; % diffusivity of species into air [m2/s]
    p.cp_liq_components = [4200 2570 4200];   % liquid specific heat [J/(kg K)]
    % we assume cp_gas_cp_N2; p.cp_gas_components = [33.46e-3;0.688e-5;0.7604e-8;-3.593e-12]./p.MW_components*1e3; % solvent (steam) specific heat in the gas phase [J/(kg K)] 
    p.eps_l_crit = [0.05 0.05 {'N/A'}];             % Drying - critical solvent content [m3_l/m3]
    p.h_M = [5e-10 5e-10 0];                   % Drying - mass transfer coefficient
    p.eps_l_eq = [0 0 1];                  % Drying - equilibrium solvent content [m3_l/m3] 
    p.latent_heat = [2260*1e3 846*1e3 {'N/A'}]; % solvent latent heat of vaporization J/kg] - model it
    %% Filtration/deliquoring/washing parameters
    p.station_diameter = 0.01;       % Filter diameter [m]
    p.Rm = 2.22e9;                   % Filter medium resistance [??m^-1];
    p.rho_sol = 1400;                % Crystal density [kg/m^3]
    p.visc_gas_phase = 1.8e-5;                % Air viscosity at room temperature and 1 atm [Pa s]
    p.surf_t = 22.39e-3;             % Surface tension [N/m]  
    p.lambda = 5;                    % Deliquoring - pore size index for deliquoring [-]
    p.k0 = 0.15;                     % Washing - diffusive backflux parameter
    p.gamma = 5;                     % Washing - diffusive back-flux parameter
    
    %% Liquid phase physical properties calculation    
    % Fixed methods
    p.mass_fr_from_conc = @(c) c./sum(c); % sum(c) is the density of the liquid phase
    p.mol_fr_from_mass_fr=@(w) w./p.MW_components./sum(w./p.MW_components);
    p.MW_liquid_phase_from_mass_fr=@(w) sum(p.MW_components.*p.mol_fr_from_mass_fr(w)); % [g/mol]
    p.rho_liquid_phase_from_mass_fr = @(w) sum(p.rho_liq_components.*w') ; % w is the mass fraction vector
    p.rho_liquid_phase_from_c=@(c) sum(c);  
    p.visc_liquid_phase_from_mass_fr= @(w) exp(sum(p.mol_fr_from_mass_fr(w').*log(p.visc_liq_components)));
    p.cp_liquid_phase_from_mass_fr = @(w) sum(p.mol_fr_from_mass_fr(w').*p.cp_liq_components);
    
    %% Calculation of cake properties 
    p.A = p.station_diameter.^2.*pi./4;                        % Filtration area [m^2]
    p.E=Porosity_function(p.station_diameter,cryst_output.CSD,cryst_output.x);           % Ouchiyama model
    p.m0=trapz(cryst_output.x,cryst_output.CSD); 
    p.m1=trapz(cryst_output.x,cryst_output.CSD.*cryst_output.x'); 
    p.m2=trapz(cryst_output.x,cryst_output.CSD.*(cryst_output.x'.^2)); 
    p.m3=trapz(cryst_output.x,cryst_output.CSD.*(cryst_output.x'.^3)); 
    p.alpha_CSD=@(dp) 180*(1-p.E)./(p.E^3*dp.^2*p.rho_sol);
    p.alpha=trapz(cryst_output.x,p.alpha_CSD(cryst_output.x).*cryst_output.CSD'/p.m0);
    p.k=1/(p.alpha*p.rho_sol*(1-p.E));      % cake permeability [?]
    p.a_V=p.m2/p.m3;
    
    %% Drying
    % Drying parameters
    % Gas phase parameters
    p.k_N2 = 0.028; % air conductivity [W/(m K)]
    p.MW_N2 = 28*1e-3; % air molecular weight [kg/mol]
    p.cp_N2 = [29e-3;0.2199e-5;0.5723e-8;-2.871e-12]/p.MW_N2*1e3; % nitrogen specific heat [J/(kg K)] 
    
    % Solid phase properties
    p.cp_s = 2200;                   % Drying - (aspirin) solid specific heat [J/(kg K)]

    % start - Antoine coefficients for water (ethanol?)
    p.A1 = 8.10765; % Antoine constant #1 % 0 to 60 C, Felder Rousseau
    p.B1 = 1750.286; % Antoine constant #2
    p.C1 = 235; % Antoine constant #3
    p.A2 = 7.96681; % 60 to 150 C
    p.B2 = 1668.21; 
    p.C2 = 228;   
    % end - Antoine coefficients for water (ethanol?)
    
    % Drying methods
    % to be adapted if you have more than one component in the liquid
    p.MWgas = @(x) (1-x(:,1))*p.MW_N2+x(:,1)*p.MW_solv; 
    %---------------------------------------------------------------
    p.cp = @(x,coeff) (coeff(1)+coeff(2)*x(:,3)+coeff(3)*x(:,3).^2+coeff(4)*x(:,3).^3); % [ J/kg K]
    p.dcpdT = @(x,coeff) (coeff(2)+2*coeff(3)*x(:,3)+3*coeff(4)*x(:,3).^2); % [ J/kg K]
    % to be adapted if you have more than one component in the liquid
    p.cp_gas = @(x) p.cp(x,p.cp_N2).*(1-x(:,1))+p.cp(x,p.cp_solventG).*x(:,1);
    p.dcp_gasdT = @(x) p.dcpdT(x,p.cp_N2).*(1-x(:,1))+p.dcpdT(x,p.cp_solventG).*x(:,1);
    p.dcp_gasdx1 =  @(x) -p.cp(x,p.cp_N2)+p.cp(x,p.cp_solventG);
%     p.dcp_gasdx2 =  @(x) p.cp(x,p.cp_N2)-p.cp(x,p.cp_solventG);
    %---------------------------------------------------------------
    p.eps_g = @(x) p.E-x(:,2);
    p.rho_g=@(x) x(:,5)./(8.314*x(:,3)).*p.MWgas(x); % gas phase density [kg/m3]
    p.h_T = @(x) 10;% ((2.19.*p.Re_bird(x).^(-2/3)).*p.rho_g(x).*p.ug(x).*p.eps_g(x).*p.cp_gas(x)./(p.cp_gas(x)*p.mu/p.k_N2).^(2/3)); % Bird, Re < 50, paragraph 14.5
    p.Psat_solv = @(x,coeff) 10.^(coeff(:,1)-coeff(:,2)./(coeff(:,3)+(x(:,3)-273.15)))*133.322; % Antoine Equation - saturation pressure [Pa]  
    
    % Mother-liquor parameters
    p.switch_antoine=@(x) 1-min(max(x(:,3)-60-273.15,0),1);
    p.coeff= @(x) [p.A1,p.B1,p.C1].*p.switch_antoine(x)+[p.A2,p.B2,p.C2].*(1-p.switch_antoine(x));   
    p.f=@(x) min((x(:,2)-p.eps_l_eq)/(p.eps_l_crit-p.eps_l_eq),1);   % Drying - drying rate efficiency function [adim] 
    p.activ_DR=@(x) min(max(x(:,4)-278,0),1);
    p.drying_driving_force = @(x) max(p.Psat_solv(x,p.coeff(x))-x(:,1)./p.MW_solv.*p.MWgas(x).*x(:,5),0);
    p.DR = @(x) p.h_M*p.a_V*p.drying_driving_force(x).*p.f(x).*p.activ_DR(x); 
    
% uncomment for calculating h_T with the correlation
%     p.Re = @(x) p.rho_g(x).*p.ug(x).*p.eps_g(x)*p.d_p/p.mu; 
%     p.Re_bird = @(x) p.Re(x)/(1-p.E);
% %     p.Sc= @(x) p.mu/(p.rho_g(x)*p.Dwat_air);
% %     p.Pr= @(x) p.cp_gas(x)*p.mu/p.k_N2;
% %     p.Nu = @(x) (2+1.3*p.Pr(x).^0.15+0.66*p.Pr(x).^0.31.*p.Re(x).^0.50);%.5; % Re is very low, so we can assume pure conduction - 0.79*(2.06/p.eps_g(x)*p.Re(x)^0.425*p.Pr(x)^(1/3));
% %     p.Nu_2 = @(x) 0.79*(2.06./p.eps_g(x).*p.Re(x).^0.425.*p.Pr(x).^(1/3));
% %     p.Sh = @(x) 1.09./p.eps_g(x).*(p.Re(x).*p.Sc(x)).^(1/3);
                 
    %     p.Vg = @(x) p.Vg_N * (101325./x(:,5)).*(x(:,3)/273.15); % flowrate of
%     drying nitrogen in [m3/s] - use it if the input is the flowrate,
%     rather than the pressure
end

