function [filtration_output,p]=model_filtration(cryst_output,p) 
    
    % Required inputs:
    % cryst_output.conc_MSMPR   = slurry concentration [kg_cryst/m3]
    % p.V_slurry                = slurry volume [mL]
    % p                         = parameters object


    %% Calculations filtration (cake formation)
    visc_liq=p.visc_liquid_phase_from_mass_fr(cryst_output.T,cryst_output.liq_mass_fr_vect); % viscosity from mass fractions, faster
%     visc_liq=exp(sum(cryst_output.liq_mass_fr_vect.*log(p.visc_liq_components(cryst_output.T)))); % viscosity from molar fractions, slower

    V_slurry_initial=p.V_slurry; %cryst_output.flowrate_MSMPR.*p.t_rot;        % Slurry volume fed at the beginning of the filtration batch [m^3]
    m_solid_initial=V_slurry_initial.*cryst_output.conc_MSMPR;    % Total solid mass filled into the filter [kg]
    V_solid_initial=m_solid_initial./p.rho_sol;        % Total solid volume filled into the filter [m^3]
    V_liq_initial=V_slurry_initial-V_solid_initial;  % Total liquid volume filled into the filter [m^3]
    V_liquid_pores_end_of_filtration=V_solid_initial*p.E/(1-p.E); % Volume of liquid in the pores at the end of filtration [m3]
    V_filt_final=V_liq_initial-V_liquid_pores_end_of_filtration;   % volume of filtrate at the end of filtration [m3]
    c=(m_solid_initial)/V_filt_final; % Mass of dry cake deposited per unit volume of filtrate [kg sol/ m3 filtrate liquid]more compact way, but same results as Brigi's

    % Determine the filtration states at the end of filtration. For having time profiles, uncomment lines 17-18                                      
    filtration_duration=visc_liq*p.alpha*c*V_filt_final^2/(2*p.A^2*p.dP)+visc_liq*p.Rm*V_filt_final/(p.A*p.dP); 
%     t_filt=0:p.time_step_filt:10*p.t_rot;
%     t_filt=t_filt(t_filt<filtration_duration);  % filtration time
    t_filt=filtration_duration;
    a_filt= visc_liq.*p.alpha.*c./(2.*p.A.^2.*p.dP);  % darcy coeff 1
    b_filt= visc_liq.*p.Rm./(p.A.*p.dP);      % darcy coeff 2
    V_filt= (-b_filt+sqrt(b_filt.^2+4.*a_filt.*t_filt))./(2.*a_filt);     % V(t) - filtrate volume [m^3]
    m_dry_cake=c*V_filt;  % mass of deposited dry cake as function of time [kg]
    H_cake=(m_dry_cake/p.rho_sol/(1-p.E))/p.A; % height of the deposited cake [m]
    dP_media_fin=p.dP./(p.alpha*m_dry_cake(end)/p.A+p.Rm)*p.Rm; % pressure drop through filter medium when cake is fully formed
    
    %% Collect outputs in the object filtration_output
    filtration_output.t_filt=t_filt;
    filtration_output.filtration_duration=filtration_duration;
    filtration_output.final_liq_mass_fr_vect=cryst_output.liq_mass_fr_vect;
    filtration_output.S_final=1;%ones(1,p.number_nodes_filt);
      
    %% Update parameters object p with cake properties calculated by this function
    p.L_cake=H_cake(end); 
    p.Pb=sum(p.pb_CSD(cryst_output.x).*cryst_output.CSD);
    p.V_liquid_pores=V_liquid_pores_end_of_filtration;
    p.m_dry_cake=m_dry_cake;
    p.dP_media_vacuum=dP_media_fin;
    rho_liq=p.rho_liquid_phase_from_mass_fr(cryst_output.liq_mass_fr_vect);
    p.S_inf=sum(0.155*(1+0.031*p.N_cap_CSD(cryst_output.x,rho_liq).^(-0.49)).*cryst_output.CSD);    % irreducible cake saturation (minimum saturation that can be achieved by displacement of the interstitial liquid by the applied vacuum    
    
end