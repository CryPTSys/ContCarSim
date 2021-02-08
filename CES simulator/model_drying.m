function drying_output = model_drying(input,drying_duration,p)

    % Required inputs:
    % input.final_liq_mass_fr_vect  =   vector liquid phase components mass fractions
    %                                   [number of components x number of nodes] or [number of components x 1] (uniform)
    % input.S_final                 =   vector of initial liquid phase saturation [1 x nodes previous model]
    % p                             =   parameters object
    %
    % Optional inputs:
    % input.nodes = vector of nodes previous model [1 x nodes previous model]
    %
    %
    % outputs are listed at the end of the code

    % Grid discretization
    p.number_nodes = p.number_nodes_drying; % local variable

    p.step_grid_drying=p.L_cake/p.number_nodes;
    p.nodes_drying=linspace(p.step_grid_drying/2,p.L_cake-p.step_grid_drying/2,...
        p.number_nodes);
    
    % Initialization + check if cake is deliquored enough for air flow
    initial_liq_mass_fr_vect=input.final_liq_mass_fr_vect;
    rho_liq=p.rho_liquid_phase_from_mass_fr(mean(initial_liq_mass_fr_vect,2));
    p.dP=p.dP_drying;
    p.S_inf=sum(0.155*(1+0.031*p.N_cap_CSD(p.x,rho_liq).^(-0.49)).*p.CSD);    % irreducible cake saturation (minimum saturation that can be achieved by displacement of the interstitial liquid by the applied vacuum    
    SR0=(mean(input.S_final)-p.S_inf)/(1-p.S_inf);   
    % if the cake is not deliquored  (SR> 20%), calculate the time needed
    % for getting SR = 20% and carry out deliquoring
    if SR0>0.2 %+0.7*p.S_inf 
        p.dP=p.dP_drying;
        p.Rm=0;
        n_switch=min(max(SR0-.2,0)*100,1);   
        thetaP0=((1-SR0)/(1.08*SR0))^(1/.88)*(1-n_switch)+((1-SR0)/(1.46*SR0))^(1/.48)*n_switch;
        residual_deliq_theta=8-thetaP0;
        visc_liq=p.visc_liquid_phase_from_mass_fr(298,mean(initial_liq_mass_fr_vect,2));
        residual_deliq=min(drying_duration,residual_deliq_theta*((p.k.*(p.dP))./(p.E.*visc_liq*(p.L_cake.^2).*(1-p.S_inf)))^-1);  % duration of deliquoring in position 5
        drying_duration=drying_duration-residual_deliq;
        input=model_deliquoring_species_grad(input,residual_deliq,p);
    end
    
    % if there is no time left after deliquoring in drying port, skip drying
    if drying_duration <= 0    
        drying_output=input;
        drying_output.final_liq_mass_fr_vect=mean(drying_output.final_liq_mass_fr_vect,2);
        drying_output.vol_cont_impurities=drying_output.final_liq_mass_fr_vect*...
            p.rho_liquid_phase_from_mass_fr(drying_output.final_liq_mass_fr_vect)./...
            p.rho_liq_components;
        drying_output.t_drying=0;
        drying_output.Tgas=298;
        drying_output.Pprofile=p.Pprofile; 
        drying_output.mass_frG=zeros(p.number_volatile_components,1);
    
    % otherwise, deliquor
    else
        
%        Initialization    
%       Retrieving the inputs and rescaling to drying grid
        initial_liq_mass_fr_vect_rescaled=zeros(p.number_components,p.number_nodes);     
        if size(initial_liq_mass_fr_vect,2)>1 && abs(initial_liq_mass_fr_vect(1,1)-initial_liq_mass_fr_vect(1,end))/initial_liq_mass_fr_vect(1)>0.01      
            for i = 1 : size(initial_liq_mass_fr_vect,1)
                initial_liq_mass_fr_vect_rescaled(i,:) = interp1(input.nodes,initial_liq_mass_fr_vect(i,:),p.nodes_drying,'pchip');
            end
        else
            initial_liq_mass_fr_vect_rescaled = repmat(initial_liq_mass_fr_vect(:,1),[1,p.number_nodes]);
        end
        initial_liq_conc=p.conc_from_mass_fr(initial_liq_mass_fr_vect_rescaled);
        S0=input.S_final;
        if size(S0,2)>1 
            S0 = interp1(input.nodes,S0,p.nodes_drying,'pchip');
        else
            S0=S0(:,1)*ones(1,p.number_nodes);
        end

        % Pre-calculations
        p.dPgdz=p.dP_drying/p.L_cake;
        p.Pprofile=linspace(101325+p.dP_drying-p.dPgdz*p.step_grid_drying/2,...
            101325-p.dPgdz*p.step_grid_drying/2,p.number_nodes);
        ug0=p.dP_drying/(p.visc_gas_phase*(p.alpha*p.rho_sol*p.L_cake*(1-p.E)+p.Rm));
        

        % Initial conditions
        epsL_0=S0*p.E;
        Tg_0=ones(1,p.number_nodes)*298;
        Ts_0=ones(1,p.number_nodes)*298;
     
        mass_frG_0=ones(1,p.number_volatile_components*p.number_nodes)*0;
        vl_0=initial_liq_conc.*repmat(epsL_0,[size(initial_liq_conc,1),1])./p.rho_liq_components;
        vL_0_vol=vl_0(1:p.number_volatile_components,:);
        vL_0_non_vol=vl_0(p.number_volatile_components+1:end,:);
        p.epsL_non_vol=sum(vL_0_non_vol,1);

        S1=triu(tril(ones(p.number_nodes)),-1); % sparsity matrix - consistent with the FVM upwind differencing scheme lower diagonal

        x0=[Ts_0';mass_frG_0';reshape(vL_0_vol',size(vL_0_vol,1)*...
            size(vL_0_vol,2),1)]';     
        S=repmat(S1,[1+2*p.number_volatile_components,1+2*p.number_volatile_components]);
        options = odeset('JPattern',S);  
        
        % solver call
        [t_drying,x]=ode15s(@drying_mex_1EB,[0 drying_duration],x0,options,...
            p.number_nodes,p.number_volatile_components,p.epsL_non_vol,...
            p.Tinlet_drying,p.cp_N2,p.rho_sol,p.E,ug0,p.Pprofile,p.MW_N2,...
            p.step_grid_drying,p.vl_crit,p.vl_eq,p.coeff_antoine',...
            p.MW_components,p.h_M,p.a_V,p.rho_liq_components,p.latent_heat,...
            p.cp_liq_components,p.cp_s,p.h_T);
        
        % output arrays preparation
        Tgas=x(:,1:p.number_nodes);
        %Tsolid=Tgas; 
        mass_frG=x(:,p.number_nodes*2+1:end-...
                (p.number_volatile_components)*p.number_nodes)';
        vol_cont_impurities=reshape(x(end,end-...
                (p.number_volatile_components)*p.number_nodes+1:end),...
                p.number_nodes,p.number_volatile_components)';
        vol_cont_impurities=[vol_cont_impurities;
            reshape(vL_0_non_vol,p.number_nodes,p.number_components-...
            p.number_volatile_components)'];
        eps_l=sum(vol_cont_impurities,1);
        final_liq_conc_vector=vol_cont_impurities./eps_l.*...
            p.rho_liq_components;
        final_liq_mass_fr_vect=final_liq_conc_vector./sum(final_liq_conc_vector);
        S_final=eps_l/p.E;
%         V_cake=pi/4*p.station_diameter^2*p.L_cake;
        rho_cake=eps_l*rho_liq+(1-p.E)*p.rho_sol;
        
        %% Collect outputs in the object drying_output
        drying_output.t_drying=t_drying;
        drying_output.vol_cont_impurities=vol_cont_impurities;
        drying_output.mass_cont_impurities=vol_cont_impurities./rho_cake.*p.rho_liq_components;
        drying_output.Tgas=Tgas;
        drying_output.Pprofile=p.Pprofile;
%         drying_output.Tsolid=Tsolid;
        drying_output.nodes=p.nodes_drying;
        drying_output.S_final=S_final;
        drying_output.final_liq_mass_fr_vect=final_liq_mass_fr_vect;
        drying_output.mass_frG=mass_frG;
        
    end 
end 