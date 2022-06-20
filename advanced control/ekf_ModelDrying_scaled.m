function x_ekf = ekf_ModelDrying_scaled(x_ekf,data,step)

    % State transition function for the drying model
    
    % Re-name variables
    u=data.u;
    Vdryer=data.Vdryer;   
    Tin=data.Tin;
    p=data.p;
    options=p.options;
    ug0=Vdryer*1e-3/60/p.A;

    % Pre-calculations
    % Pressure filtration
    p.dPgdz=(u.P_compr-8700)/p.L_cake;
    p.Pprofile=linspace(101325+u.P_compr-p.dPgdz*p.step_grid_drying/2,...
        101325-p.dPgdz*p.step_grid_drying/2,p.number_nodes);
    
    % Initial conditions - Troom=22°C
    Tg_0=x_ekf(1:p.number_nodes)*300;
    Ts_0=x_ekf(1+p.number_nodes:2*p.number_nodes)*300;
    mass_frG_0=x_ekf(2*p.number_nodes+1:3*p.number_nodes)*0.01;    
    vl_0=x_ekf(3*p.number_nodes+1:4*p.number_nodes)*0.01;
    p.epsL_non_vol=sum(vl_0);
    
    x0=[Tg_0; Ts_0; mass_frG_0; vl_0]';   
    
    cp_g=p.cp(Tin);
    rho_g=101325/8.314/(Tin)*0.02897;
    
    p.h_M=ug0*rho_g/p.zeta/47.2e-6/p.a_V/1e5;
    p.h_T=cp_g*rho_g*ug0/p.a_V/47.2e-6/p.zeta;  %(1+randn(1)*0.01);

    % solver call
    [~,x]=ode15s(@drying_model,[0 step],x0,options,...
        p.number_nodes,1,zeros(1,p.number_nodes),...
        Tin,p.cp_gas_components',p.rho_sol,p.E,ug0,p.Pprofile,p.MW_air,...
        p.step_grid_drying,p.vl_crit,p.vl_eq,p.coeff_antoine',...
        p.MW_components,p.h_M,p.a_V,p.rho_liq_components,p.latent_heat,...
        p.cp_liq_components,p.cp_s,p.h_T);

    x_ekf=x(end,:)';
    x_ekf=x_ekf./([300*ones(2*p.number_nodes,1); 0.01*ones(2*p.number_nodes,1)]);

end 