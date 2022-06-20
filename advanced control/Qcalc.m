%% Calculation of process noise covariance

function Q = Qcalc(x, data) % consider supplying filter and x_Estim: when close to equilibrium, stop updating Q
    VarCovMatrix =  zeros(4);
    
    VarCovMatrix(1,1)=0.002^2;
    VarCovMatrix(4,4)=(7e-3)^2;
    VarCovMatrix(2,2)=(4e-11)^2;
    VarCovMatrix(3,3)=0.005^2;

    Tin=data.Tin;
    p=data.p;       
    cp_g=p.cp(Tin);
    rho_g=101325/8.314/(Tin)*0.02897;
    Vdryer=data.Vdryer;
    
    ug0=Vdryer*1e-3/60/p.A;
    pars0(1)=ug0;
    pars0(2)=ug0*rho_g/p.zeta/47.2e-6/p.a_V/1e5;
    pars0(3)=cp_g*rho_g*ug0/p.a_V/47.2e-6/p.zeta;
    pars0(4)=0.35;

    J=jacobianest(@fun,pars0);
    Q=J*VarCovMatrix*J';

    function dxdt = fun(pars0)
        % State transition function for the drying model         
        u=data.u;
        Tin=data.Tin;
        p=data.p;

        % uncertain parameters
        ug0=pars0(1);
        p.h_M=pars0(2);
        p.h_T=pars0(3);
        p.E=pars0(4);

        % Pre-calculations
        % Pressure filtration
        p.dPgdz=(u.P_compr-8700)/p.L_cake;
        p.Pprofile=linspace(101325+u.P_compr-p.dPgdz*p.step_grid_drying/2,...
            101325-p.dPgdz*p.step_grid_drying/2,p.number_nodes);

        % Initial conditions - Troom=22°C
        Tg_0=x(1:p.number_nodes);
        Ts_0=x(1+p.number_nodes:2*p.number_nodes);
        mass_frG_0=x(2*p.number_nodes+1:3*p.number_nodes);    
        vl_0=x(3*p.number_nodes+1:4*p.number_nodes);
        p.epsL_non_vol=sum(vl_0);

        x0=[Tg_0'*300; Ts_0'*300; mass_frG_0'*0.01; vl_0'*0.01]';   

        % solver call
        dxdt=drying_model(1,x0,...
            p.number_nodes,1,zeros(1,p.number_nodes),...
            Tin,p.cp_gas_components',p.rho_sol,p.E,ug0,p.Pprofile,p.MW_air,...
            p.step_grid_drying,p.vl_crit,p.vl_eq,p.coeff_antoine',...
            p.MW_components,p.h_M,p.a_V,p.rho_liq_components,p.latent_heat,...
            p.cp_liq_components,p.cp_s,p.h_T);

        dxdt=dxdt./([300*ones(2*p.number_nodes,1); 0.01*ones(2*p.number_nodes,1)]);

    end
end