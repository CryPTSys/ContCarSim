function d = model_deliquoring_design_charts(p)

    %% Deliquoring model calculation - solution obtained through design charts
    t_deliq=0:.1:p.t_deliq_final;
    ThetaP=(t_deliq.*p.k.*p.dP)./(p.E.*p.visc.*(p.L_cake.^2).*(1-p.S_inf)); % Dimensionless deliquoring time [-]
    n_switch=min(max(ThetaP-1.915,0)*10,1);           
    SR=1./(1+1.46.*ThetaP.^0.48).*n_switch+1./(1+1.08.*ThetaP.^0.88).*(1-n_switch);
    S=p.S_inf+SR.*(1-p.S_inf);   
%     V_liq_pores_deliq=S*filt_output.V_liquid_pores_filtration; % Volume of liquid in the pores during deliq as function of time
%     V_deliq=(1-S)*filt_output.V_liquid_pores_filtration;   % Volume of filtrate during deliq as function of time
%     V_liq_pores_eq=p.S_inf*filt_output.V_liquid_pores_filtration;  % Volume of liquid in the pores at eq   
%     Volume_cake=filt_output.V_liquid_pores_filtration+filt_output.m_dry_cake(end)/p.rho_sol;
    solvent_content_vol_deliq=S*p.E;%V_liq_pores_deliq/Volume_cake; % Volumetric fraction of solvent content during deliq
    solvent_content_vol_eq=p.S_inf*p.E;%V_liq_pores_eq/Volume_cake; % Volumetric fraction of solvent content at mech eq

    %% Collect outputs in the object d
    d.t_deliq=t_deliq;
%     d.V_deliq=V_deliq;
    d.solvent_content_vol_eq=solvent_content_vol_eq;
    d.solvent_content_vol_deliq=solvent_content_vol_deliq;
    d.S=S;
    d.p.S_inf=p.S_inf;
    d.ThetaP=ThetaP;
    d.SR=SR;
    
end