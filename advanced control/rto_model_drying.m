function drying_output = rto_model_drying(input,drying_duration,p)

    % Required inputs:
    % input.S_final                 =   initial liquid phase saturation [average]
    % p                             =   parameters object
    % outputs are listed at the end of the code

    % Grid discretization
    number_nodes = 10;
    step_grid_drying=p.L_cake/number_nodes;

    % Initialization + check if cake is deliquored enough for air flow
    p.S_inf=0.085;   
    SR0=(mean(input.S_final)-p.S_inf)/(1-p.S_inf);   
    % if the cake is not deliquored  (SR> 20%), calculate the time needed
    % for getting SR = 20% and carry out deliquoring
    if SR0>0.2 
        n_switch=min(max(SR0-.2,0)*100,1);   
        thetaP0=((1-SR0)/(1.08*SR0))^(1/.88)*(1-n_switch)+((1-SR0)/(1.46*SR0))^(1/.48)*n_switch;
        residual_deliq_theta=8-thetaP0;
        visc_liq=p.visc_liq_components(p.T_room);
        residual_deliq=min(drying_duration,residual_deliq_theta*...
            ((p.k.*(p.dP))./(p.E.*visc_liq*(p.L_cake.^2).*(1-p.S_inf)))^-1);  
        drying_duration=drying_duration-residual_deliq;
        input=rto_model_deliquoring(input,residual_deliq,p);
    end
    
    % if there is no time left after deliquoring in drying port, skip drying
    if drying_duration <= 0    
        drying_output=input;
        drying_output.vol_cont_impurities=drying_output.final_liq_mass_fr_vect*...
            p.rho_liquid_phase_from_mass_fr(drying_output.final_liq_mass_fr_vect)./...
            p.rho_liq_components;  
    % otherwise, dry
    else
        
%       Initialization    
        initial_liq_conc=842;
        S0=input.S_final;
        S0=mean(S0)*ones(1,p.number_nodes);

        
        % Pressure filtration
        dPgdz=(p.dP-p.dP_media)/p.L_cake; %
        Pprofile=linspace(101325+p.dP-dPgdz*step_grid_drying/2,...
            101325-dPgdz*step_grid_drying/2,number_nodes);
        ug0=p.dP/(p.visc_gas_phase*(p.alpha*p.rho_sol*p.L_cake*(1-p.E)+p.Rm));
        
        % Initial conditions - Troom=22°C
        epsL_0=S0*p.E;
        Tg_0=ones(1,p.number_nodes)*295.25;
        Ts_0=ones(1,p.number_nodes)*295.25; 
        Tin_C=Ts_0(1)-273.15;
        
        mass_frG_0=ones(1,p.number_volatile_components*p.number_nodes)*0;
        vl_0=initial_liq_conc.*repmat(epsL_0,[size(initial_liq_conc,1),1])./p.rho_liq_components;
        vL_0_vol=vl_0(1:p.number_volatile_components,:);
        vL_0_non_vol=vl_0(p.number_volatile_components+1:end,:);
        epsL_non_vol=sum(vL_0_non_vol,1);
        
        Tgas(1,:)=Tg_0;
        Tsolid(1,:)=Tg_0;
        vol_cont_impurities(1,:)=vl_0 ;
        eps_l(1,:)=vol_cont_impurities(end,:);
        rho_cake=(eps_l(1,:)*p.rho_liq_components+(1-p.E)*p.rho_sol);
        mass_avg(1)=mean(vol_cont_impurities(1,:)./rho_cake.*p.rho_liq_components);  
        t_drying(1)=0;
        mass_frG(1,:)=mass_frG_0;       

        S1=triu(tril(ones(p.number_nodes)),-1); % sparsity matrix - consistent with the FVM upwind differencing scheme lower diagonal

        x0=[Tg_0'; Ts_0'; mass_frG_0'; vL_0_vol']';%reshape(vL_0_vol',size(vL_0_vol,1)*...
%             size(vL_0_vol,2),1)]';     
        S=repmat(S1,[2+2*p.number_volatile_components,2+2*p.number_volatile_components]);
        options = odeset('JPattern',S);  
        
        % inlet temperature profile calculation
        heat_loss_pars=[0.699889288242627,68.884855712110380,-0.71788879977705];%[0.67,130,-0.35]; % 1 2 3 7 8 9 |33.19
        a_hl=heat_loss_pars(1);
        tau=heat_loss_pars(2);
        ug_exp=heat_loss_pars(3);
        a=.155./(tau+.155);       
        time_steps_HL=linspace(0,drying_duration,max(ceil(drying_duration/.155),2));
        time_steps_simulation=linspace(0,drying_duration,max(ceil(drying_duration/10),2));
        for i =1:length(time_steps_HL)
            Tin_C=(1-a)*Tin_C+a*(p.Tinlet_drying-273.15-(a_hl+ug0*...
              ug_exp)*(p.Tinlet_drying-295));              
            Tin(i)=Tin_C+273.15;
        end
        Tin=interp1(time_steps_HL,Tin, time_steps_simulation);
        
        % simulation
      
        for i = 1:1:length(time_steps_simulation)-1       
            
            % calculation of hT and mT coefficients
            cp_g=p.cp(Tin(i));
            rho_g=101325/8.314/(Tin(i))*0.02897;            
            p.h_M=ug0*rho_g/p.zeta/47.2e-6/p.a_V/1e5;
            p.h_T=cp_g*rho_g*ug0/p.a_V/47.2e-6/p.zeta;
            
            % solver call
            [tt,x]=ode15s(@drying_model,[time_steps_simulation(i) time_steps_simulation(i+1)],x0,options,...
                number_nodes,p.number_volatile_components,epsL_non_vol,...
                Tin(i),p.cp_gas_components',p.rho_sol,p.E,ug0,Pprofile,p.MW_air,...
                step_grid_drying,p.wl_crit,p.wl_eq,p.coeff_antoine',...
                p.MW_components,p.h_M,p.a_V,p.rho_liq_components,p.latent_heat,...
                p.cp_liq_components,p.cp_s,p.h_T);

            % output arrays
            x0=x(end,:);
            vol_cont_impurities(i+1,1:p.number_nodes*p.number_volatile_components)=...
                x(end,end-(p.number_volatile_components)*p.number_nodes+1:end);
            eps_l(i+1,:)=vol_cont_impurities(end,:);
            rho_cake=(eps_l(i+1,:)*p.rho_liq_components+(1-p.E)*p.rho_sol);
            mass_avg(i+1)=mean(vol_cont_impurities(i+1,:)./rho_cake.*p.rho_liq_components);      
            vol_cont_impurities(i+1,p.number_nodes*p.number_volatile_components+1:end)=...
                vL_0_non_vol;  
        end

        
        %% Collect outputs in the object drying_output
        drying_output.mass_avg_profile=mass_avg;
        
    end 
end 