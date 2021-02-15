function [x,y] = model_drying(cycle_time,Dt,p,u,x,y,n_batch,pos)
    % differential drying model
    
    %%  Inputs list
    %   cycle_time = cycle timer
    %   Dt = duration of deliquoring step
    %   p = properties object
    %   x = states (+additional properties) object
    %   y = measurements vector
    %   n_cycle = cycle number
    %   pos = port in which filtration is occurring               

    %% update deliq. equilibrium saturation to current pressure drop
    rho_liq=p.rho_liquid_phase_from_mass_fr(mean(x.pos4.liq_mass_fr_vect,2));
    x.(['pos' num2str(pos)]).S_inf=sum(0.155*(1+...
        0.031*p.N_cap_CSD(x.(['pos' num2str(pos)]).x_perc,rho_liq,...
        x.(['pos' num2str(pos)]).E,u.dP_drying,x.(['pos' num2str(pos)]).L_cake).^...
        (-0.49)).*x.(['pos' num2str(pos)]).CSD_perc);   
    
    %% Initial conditions (creation of local objects with shorter names for conciseness)
    t=cycle_time;
    xx=x.(['pos' num2str(pos)]);        
    initial_liq_mass_fr_vect=xx.liq_mass_fr_vect;    
    SR0=(mean(xx.S)-xx.S_inf)/(1-xx.S_inf);
    residual_deliq=0;
    
    %% Check if the cake is not deliquored enough for air to flow (SR> 20%)
    % if that is the case, estimate the time needed
    % for getting SR = 20% using design charts and simulate deliquoring with
    % PDE model
    if SR0>0.2 % deliquoring happens
        n_switch=min(max(SR0-.2,0)*100,1);   
        thetaP0=((1-SR0)/(1.08*SR0))^(1/.88)*(1-n_switch)+((1-SR0)/(1.46*SR0))^(1/.48)*n_switch;
        residual_deliq_theta=8-thetaP0; % estimation of deliquoring time for reaching SR=20% using design charts
        visc_liq=p.visc_liquid_phase_from_mass_fr(xx.T,mean(initial_liq_mass_fr_vect,2));
        residual_deliq=min(Dt,residual_deliq_theta*((xx.k.*(u.dP_drying))./(xx.E.*visc_liq*(xx.L_cake.^2).*(1-xx.S_inf)))^-1);  % duration of deliquoring in position 4
        Dt=Dt-residual_deliq;
        u.dP=u.dP_drying;
        [x,y]=model_deliquoring_species_grad(t,residual_deliq,p,u,x,y,n_batch,pos); % deliquoring simulation
        % measurement vectors update
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_batch)]).t_drying=...
            y.(['pos' num2str(pos)]).(['cycle_' num2str(n_batch)]).t_filt;
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_batch)]).wg=...
            zeros(p.number_volatile_components,length(y.(['pos' num2str(pos)]).(['cycle_' num2str(n_batch)]).t_filt));
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_batch)]).Tg=298*...
            ones(1,length(y.(['pos' num2str(pos)]).(['cycle_' num2str(n_batch)]).t_filt));
    end
    
    %% Drying simulation
    if Dt > 0    % drying occurs
        %Initial conditions (creation of local objects with shorter names for conciseness)
        xx=x.(['pos' num2str(pos)]);
        initial_liq_mass_fr_vect=xx.liq_mass_fr_vect;       
        S0=xx.S;
        
        % Rescale initial saturation and composition profiles to match current grid
        if size(initial_liq_mass_fr_vect,2)>1 && abs(initial_liq_mass_fr_vect(1,1)-initial_liq_mass_fr_vect(1,end))/initial_liq_mass_fr_vect(1)>1e-6     
            initial_liq_mass_fr_vect_rescaled=zeros(p.number_components,xx.number_nodes_drying);
            for i = 1 : size(initial_liq_mass_fr_vect,1)
                initial_liq_mass_fr_vect_rescaled(i,:) = interp1(xx.nodes,initial_liq_mass_fr_vect(i,:),xx.nodes_drying,'pchip');
            end
        else
            initial_liq_mass_fr_vect_rescaled = ones(size(initial_liq_mass_fr_vect,1),...
                xx.number_nodes_drying).*initial_liq_mass_fr_vect(:,1);
        end
        initial_liq_conc=p.conc_from_mass_fr(initial_liq_mass_fr_vect_rescaled);
        if size(S0,2)>1 && abs(S0(1,1)-S0(1,end))/S0(1)>1e-5 
            S0 = interp1(xx.nodes,S0,xx.nodes_drying);
        else
            S0 = ones(1,xx.number_nodes_drying).*S0(:,1);
        end

        % Calculations
        dPgdz=(u.dP_drying-xx.dP_mesh)/xx.L_cake;        
        Pgin=101325+u.dP_drying-dPgdz*xx.step_grid_drying/2; % pressure first node
        Pgout=101325-dPgdz*xx.step_grid_drying/2;          % pressure last node
        Pprofile=linspace(Pgin,Pgout,xx.number_nodes_drying);
        ug0=u.dP_drying/(p.visc_gas_phase*(xx.alpha*p.rho_sol*xx.L_cake*(1-xx.E)+p.Rm)); % temperature and saturation dependencies added in mex file        
        epsL_0=S0*xx.E; % liquid phase volumetric fraction in cake
        Tg_0=xx.Tg; % gas temperature profile
        mass_frG_0=xx.gas_mass_fr_vect;    % gas phase composition - mass fractions
        vl_0=initial_liq_conc.*repmat(epsL_0,[size(initial_liq_conc,1),1])./p.rho_liq_components; % liquid phase composition - volumetric fractions
        vL_0_vol=vl_0(1:p.number_volatile_components,:); % liquid phase composition - volumetric fractions - only volatile componentes
        vL_0_non_vol=vl_0(p.number_volatile_components+1:end,:); % liquid phase composition - volumetric fractions - only non-volatile componentes
        epsL_non_vol=sum(vL_0_non_vol,1); % volumetric fraction of non-volatile lliquid components in cake
        x0=[Tg_0'; mass_frG_0'; reshape(vL_0_vol',size(vL_0_vol,1)*...
            size(vL_0_vol,2),1)]'; % initial state
        
        % Time vector
        t_drying=unique([0 (p.drying_sampling_time-round(rem(t,p.drying_sampling_time),6)):p.drying_sampling_time:Dt Dt]);   
                
        % sparsity matrix - consistent with the FVM upwind differencing scheme
        S1=triu(tril(ones(xx.number_nodes_drying)),-1); % lower diagonal
        S=repmat(S1,[1+2*p.number_volatile_components,1+2*p.number_volatile_components]);
 
        options = odeset('JPattern',S);
        
        % solver call
        [t_drying,output]=ode15s(@drying_mex_1EB,t_drying,x0,options,...
            xx.number_nodes_drying,p.number_volatile_components,epsL_non_vol,...
            u.Tinlet_drying,p.cp_air,p.rho_sol,xx.E,ug0,Pprofile,p.MW_air,...
            xx.step_grid_drying,p.vl_crit,p.vl_eq,p.coeff_antoine',...
            p.MW_components,p.h_M,xx.a_V,p.rho_liq_components,p.latent_heat,...
            p.cp_liq_components,p.cp_s, p.h_T);
        
        % if t_drying has only two elements, ode15s gives output with more time steps
        % the next if loop is necessary for having the output consistent with t_drying    
        if length(t_drying)==2
            output=[output(1,:); output(end,:)];
        end
        
        % Extract states from solver output
        Tgas=output(:,1:xx.number_nodes_drying);
        Tsolid=Tgas;
        mass_frG=output(:,xx.number_nodes_drying+1:end-...
                (p.number_volatile_components)*xx.number_nodes_drying);
        vol_cont_impurities=max(reshape(output(end,end-...
                (p.number_volatile_components)*xx.number_nodes_drying+1:end),...
                xx.number_nodes_drying,p.number_volatile_components)',p.vl_eq(1:p.number_volatile_components));
        vol_cont_impurities=[vol_cont_impurities;
            reshape(vL_0_non_vol,xx.number_nodes_drying,p.number_components-...
            p.number_volatile_components)'];
        eps_l=sum(vol_cont_impurities,1);
        final_liq_conc_vector=vol_cont_impurities./eps_l.*...
            p.rho_liq_components;
        final_liq_mass_fr_vect=final_liq_conc_vector./sum(final_liq_conc_vector);
        S_final=eps_l/xx.E;
        
        %% Updates states vector (x) and measurement vector (y)
        % states      
        t_drying=t_drying+t+residual_deliq; % switch to process time
        x.(['pos' num2str(pos)]).Tg=Tgas(end,:);
        x.(['pos' num2str(pos)]).Ts=Tsolid(end,:);
        x.(['pos' num2str(pos)]).S=S_final(end,:);
        x.(['pos' num2str(pos)]).nodes=xx.nodes_deliq;
        x.(['pos' num2str(pos)]).liq_mass_fr_vect=final_liq_mass_fr_vect;
        x.(['pos' num2str(pos)]).gas_mass_fr_vect=mass_frG(end,:);

        % measurements
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_batch)]).t_drying=...
            [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_batch)]).t_drying,...
            t_drying(2:end)'];
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_batch)]).Tg=...
            [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_batch)]).Tg, Tgas(2:end,end)'];
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_batch)]).wg=...
            [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_batch)]).wg,...
            mass_frG(2:end,[1:p.number_volatile_components]*xx.number_nodes_drying)'];
    end
end
    
 