function d = model_deliquoring_design_charts_eqCSD(cryst_output,filt_output,p)

    %% Deliquoring model calculation - solution obtained through design charts
    t_deliq=0:.1:p.t_deliq_final;
    k_av= 1/(filt_output.alpha*p.rho_sol*(1-filt_output.E));      % cake permeability [?]
    L_cake=filt_output.H_cake(end);                   % cake height at the end of filtration [m]
    %N_cap=(filt_output.E^3*(filt_output.m1/filt_output.m0)^2*(p.rho_liq*9.81.*L_cake+p.dP))/((1-filt_output.E)^2*L_cake.*p.surf_t);
    N_cap_CSD=@(x) (filt_output.E^3*(x).^2*(p.rho_liq*9.81.*L_cake+p.dP))/((1-filt_output.E)^2*L_cake.*p.surf_t);
    N_cap=trapz(cryst_output.x,N_cap_CSD(cryst_output.x).*cryst_output.CSD'/filt_output.m0);
    S_inf= 0.155.*(1+0.031.*N_cap.^(-0.49));   % irreducible cake saturation (minimum saturation that can be achieved by displacement of the interstitial liquid by the applied vacuum
    ThetaP=(t_deliq.*k_av.*p.dP)./(filt_output.E.*p.visc.*(L_cake.^2).*(1-S_inf)); % Dimensionless deliquoring time [-]
    n_switch=min(max(ThetaP-1.915,0)*10,1);           
    SR=1./(1+1.46.*ThetaP.^0.48).*n_switch+1./(1+1.08.*ThetaP.^0.88).*(1-n_switch);
    S=S_inf+SR.*(1-S_inf);   
    V_liq_pores_deliq=S*filt_output.V_liquid_pores_filtration; % Volume of liquid in the pores during deliq as function of time
    V_deliq=(1-S)*filt_output.V_liquid_pores_filtration;   % Volume of filtrate during deliq as function of time
    V_liq_pores_eq=S_inf*filt_output.V_liquid_pores_filtration;  % Volume of liquid in the pores at eq   
    Volume_cake=filt_output.V_liquid_pores_filtration+filt_output.m_dry_cake(end)/p.rho_sol;
    solvent_content_vol_deliq=V_liq_pores_deliq/Volume_cake; % Volumetric fraction of solvent content during deliq
    solvent_content_vol_eq=V_liq_pores_eq/Volume_cake; % Volumetric fraction of solvent content at mech eq

    %% Collect outputs in the object d
    d.t_deliq=t_deliq;
    d.V_deliq=V_deliq;
    d.solvent_content_vol_eq=solvent_content_vol_eq;
    d.solvent_content_vol_deliq=solvent_content_vol_deliq;
    d.S=S;
    d.S_inf=S_inf;
    d.ThetaP=ThetaP;
    d.SR=SR;
    
end