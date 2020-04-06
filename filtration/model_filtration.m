function f=model_filtration(t,cryst_output,p) 

    %% Calculations filtration (cake formation)
    A = p.Filter_d.^2.*pi./4;                        % Filtration area [m^2]
    E=Porosity_function(p.Filter_d,cryst_output.CSD,cryst_output.x);           % Ouchiyama model
    eps_s= 1-E;                                    % Solid volume fraction in cake [m^3 cryst/m^3 cake]
    m0=trapz(cryst_output.x,cryst_output.CSD); 
    m1=trapz(cryst_output.x,cryst_output.CSD.*cryst_output.x'); 

    V_slurry_initial=cryst_output.flowrate_MSMPR.*p.t_rot;        % Slurry volume fed at the beginning of the filtration batch [m^3]
    %L_slurry_init=V_slurry_in./A;                 % Height of the slurry column at the beginning [m]

    m_solid_initial=V_slurry_initial.*cryst_output.conc_MSMPR;    % Total solid mass filled into the filter [kg]
    V_solid_initial=m_solid_initial./p.rho_sol;        % Total solid volume filled into the filter [m^3]
    w_solid_initial=cryst_output.conc_MSMPR./(cryst_output.conc_MSMPR+p.rho_liq.*(1-cryst_output.conc_MSMPR./p.rho_sol));  % Solid mass fraction in the initial slurry [-]

    mass_slurry_initial=m_solid_initial/w_solid_initial;            % Total cake mass filled into the filter [kg]

    V_liq_initial=V_slurry_initial-V_solid_initial;  % Total liquid volume filled into the filter [m^3]

    V_liquid_pores_filtration=V_solid_initial*E/(1-E); % at the end of filtration [m3]
    V_filt_final=V_liq_initial-V_liquid_pores_filtration;   % volume of filtrate at the end of filtration [m3]

    % c_brig= 1./(((1-w_solid_initial)./(w_solid_initial.*p.rho_liq))-...
    %     ((E)./(eps_s.*p.rho_sol)));   % Mass of dry cake deposited per unit volume of filtrate [kg cryst/ m3 filtrate liquid]

    c=(m_solid_initial)/V_filt_final; % Mass of dry cake deposited per unit volume of filtrate [kg sol/ m3 filtrate liquid]more compact way, but same results as Brigi's

    % Cake resistance
    alpha_CSD=@(dp) 180*(1-E)./(E^3*dp.^2*p.rho_sol);
    alpha=trapz(cryst_output.x,alpha_CSD(cryst_output.x).*cryst_output.CSD'/m0);

    % Determine the filtration states as function of time                                       
    t_filt_total=p.visc*alpha*c*V_filt_final^2/(2*A^2*p.dP)+p.visc*p.Rm*V_filt_final/(A*p.dP); 
    t_filt=t(t<t_filt_total);  % filtration time
    a_filt= p.visc.*alpha.*c./(2.*A.^2.*p.dP);  % darcy coeff 1
    b_filt= p.visc.*p.Rm./(A.*p.dP);      % darcy coeff 2
    V_filt= (-b_filt+sqrt(b_filt.^2+4.*a_filt.*t_filt))./(2.*a_filt);     % V(t) - filtrate volume [m^3]
    m_dry_cake=c*V_filt;  % mass of deposited dry cake as function of time [kg]
    m_wet_cake=m_dry_cake*(1+E/(1-E)/p.rho_sol*p.rho_liq);
    Q=(alpha*p.visc*c*V_filt/(p.dP*A^2)+p.Rm*p.visc/(p.dP*A)).^-1;
    H_cake=(m_dry_cake/p.rho_sol/(1-E))/A; % height of the deposited cake [m]
    m_liq_pores=m_wet_cake-m_dry_cake; % mass of liquid in the pores during filtration [m3]
    solvent_content_mass=m_liq_pores./(m_liq_pores+m_dry_cake);
    solvent_content_vol_filt=m_liq_pores/p.rho_liq./(m_liq_pores/p.rho_liq+m_dry_cake/p.rho_sol);

    %% Collect outputs in the object f
    f.solvent_content_mass=solvent_content_mass;
    f.solvent_content_vol_filt=solvent_content_vol_filt;
    f.t_filt=t_filt;
    f.V_filt=V_filt;
    f.t_filt_total=t_filt_total;
    f.H_cake=H_cake;
    f.alpha=alpha;
    f.E=E;
    f.m0=m0;
    f.m1=m1;
    f.V_liquid_pores_filtration=V_liquid_pores_filtration;
    f.m_dry_cake=m_dry_cake;
end