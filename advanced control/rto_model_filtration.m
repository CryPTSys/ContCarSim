function [filtration_output,p]=rto_model_filtration(cryst_output,p) 
    
    %% Calculations filtration (cake formation)
    visc_liq=p.visc_liq_components(p.T_room); % viscosity from mass fractions, faster
    V_slurry_initial=p.V_slurry; %cryst_output.flowrate_MSMPR.*p.t_rot;        % Slurry volume fed at the beginning of the filtration batch [m^3]
    m_solid_initial=V_slurry_initial.*cryst_output.conc_slurry;    % Total solid mass filled into the filter [kg]
    V_solid_initial=m_solid_initial./p.rho_sol;        % Total solid volume filled into the filter [m^3]
    V_liq_initial=V_slurry_initial-V_solid_initial;  % Total liquid volume filled into the filter [m^3]
    V_liquid_pores_end_of_filtration=V_solid_initial*p.E/(1-p.E); % Volume of liquid in the pores at the end of filtration [m3]
    V_filt_final=V_liq_initial-V_liquid_pores_end_of_filtration;   % volume of filtrate at the end of filtration [m3]
    c=(m_solid_initial)/V_filt_final; % Mass of dry cake deposited per unit volume of filtrate [kg sol/ m3 filtrate liquid]more compact way, but same results as Brigi's
 
    % Determine the filtration states at the end of filtration. For having time profiles, uncomment lines 17-18                                      
    filtration_duration=visc_liq*p.alpha*c*V_filt_final^2/(2*p.A^2*p.dP)+visc_liq*p.Rm*V_filt_final/(p.A*p.dP); 
    t_filt=filtration_duration;
    a_filt= visc_liq.*p.alpha.*c./(2.*p.A.^2.*p.dP);  % darcy coeff 1
    b_filt= visc_liq.*p.Rm./(p.A.*p.dP);      % darcy coeff 2
    V_filt= (-b_filt+sqrt(b_filt.^2+4.*a_filt.*t_filt))./(2.*a_filt);     % V(t) - filtrate volume [m^3]
    m_dry_cake=c*V_filt;  % mass of deposited dry cake as function of time [kg]    
    dP_media_fin=p.dP./(p.alpha*m_dry_cake(end)/p.A+p.Rm)*p.Rm; % pressure drop through filter medium when cake is fully formed
    
    %% Collect outputs in the object filtration_output
    filtration_output.filtration_duration=filtration_duration;
    filtration_output.S_final=1;    
      
    %% Update parameters object p with cake properties calculated by this function
    p.L_cake=p.V_slurry*cryst_output.conc_slurry/p.rho_sol/(1-p.E)/p.A;
    if p.L_cake < 0
        disp(1)
    end
    p.Pb=4.3149e3; 
    p.dP_media=dP_media_fin;
    
end