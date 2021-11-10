function [x,y] = model_drying(batch_time,Dt,p,d,u,x,y,n_batch,pos)
    %% differential drying model    
    %  Inputs list
    %   batch_time = cycle timer
    %   Dt = duration of deliquoring step
    %   p = properties object
    %   x = states (+additional properties) object
    %   y = measurements vector
    %   n_cycle = cycle number
    %   pos = station in which filtration is occurring               

    %% update deliq. equilibrium saturation to current pressure drop
    rho_liq=p.rho_liq_components;
    x.(['pos' num2str(pos)]).S_inf=0.085;
    
    %% Initial conditions (creation of local objects with shorter names for conciseness)
    t=batch_time;
    xx=x.(['pos' num2str(pos)]);        
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
        visc_liq=p.visc_liq_components(xx.T);
        residual_deliq=min(Dt,residual_deliq_theta*((xx.k.*(u.dP))./(xx.E.*visc_liq*(xx.L_cake.^2).*(1-xx.S_inf)))^-1);  % duration of deliquoring in position 4
        Dt=Dt-residual_deliq;
        u.dP=u.dP;
        [x,y]=model_deliquoring_grad(t,residual_deliq,p,u,x,y,n_batch,pos); % deliquoring simulation
        % measurement vectors 
        y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).w_EtOH_gas=...
            zeros(p.number_volatile_components,length(y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).t));
        y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).Tg_top=...
            y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).Tg_top(end)*...
            ones(1,length(y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).t));
        y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).Vdryer=...
            zeros(1,length(y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).t));         
        y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).Tg_bot=...
            y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).Tg_bot(end)*...
            ones(1,length(y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).t));
        y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).Ts_bot=...
            y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).Ts_bot(end)*...
            ones(1,length(y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).t));    
    end
    
    %% Drying simulation
    if Dt > 0    % drying occurs
        %Initial conditions (creation of local objects with shorter names for conciseness)
        xx=x.(['pos' num2str(pos)]);      
        S0=xx.S;
        
        % Rescale initial saturation and composition profiles to match current grid
        if size(S0,2)>1 && abs(S0(1,1)-S0(1,end))/S0(1)>1e-5 
            S0 = interp1(xx.nodes,S0,xx.nodes_drying);
        else
            S0 = ones(1,xx.number_nodes_drying).*S0(:,1);
        end                
        initial_liq_conc=p.rho_liq_components;

        % Calculations
        dPgdz=(u.dP-xx.dP_mesh)/xx.L_cake;        
        Pgin=101325+u.dP-dPgdz*xx.step_grid_drying/2; % pressure first node
        Pgout=101325-dPgdz*xx.step_grid_drying/2;          % pressure last node
        Pprofile=linspace(Pgin,Pgout,xx.number_nodes_drying);
        ug0=u.dP/(p.visc_gas_phase*(xx.alpha*p.rho_sol*xx.L_cake*(1-xx.E)+p.Rm(pos))); % temperature and saturation dependencies added in mex file        
        epsL_0=S0*xx.E; % liquid phase volumetric fraction in cake
        Tg_0=xx.Tg; % gas temperature profile
        Ts_0=xx.Ts;
        mass_frG_0=xx.gas_mass_fr_vect;    % gas phase composition - mass fractions
        vl_0=initial_liq_conc.*repmat(epsL_0,[size(initial_liq_conc,1),1])./p.rho_liq_components; % liquid phase composition - volumetric fractions
        vL_0_vol=vl_0(1:p.number_volatile_components,:); % liquid phase composition - volumetric fractions - only volatile componentes
        vL_0_non_vol=vl_0(p.number_volatile_components+1:end,:); % liquid phase composition - volumetric fractions - only non-volatile componentes
        epsL_non_vol=sum(vL_0_non_vol,1); % volumetric fraction of non-volatile lliquid components in cake
        x0=[Tg_0'; Ts_0'; mass_frG_0'; reshape(vL_0_vol',size(vL_0_vol,1)*...
            size(vL_0_vol,2),1)]'; % initial state
        
        % Time vector
        t_drying=unique([0 (p.drying_sampling_time-round(rem(t,p.drying_sampling_time),6)):p.drying_sampling_time:Dt Dt]);   
                
        % sparsity matrix - consistent with the FVM upwind differencing scheme
        S1=triu(tril(ones(xx.number_nodes_drying)),-1); % lower diagonal
        S=repmat(S1,[2+2*p.number_volatile_components,2+2*p.number_volatile_components]);
 
        options = odeset('JPattern',S);
        
        % update inlet temperature after HL
        heat_loss_pars=[0.6998,68.88485571,-0.7178];   
        a_hl=heat_loss_pars(1);
        tau=heat_loss_pars(2);
        ug_exp=heat_loss_pars(3);
        a=.155./(tau+.155);               
        time_steps_HL=linspace(0,t_drying(end),max(ceil(t_drying(end)/.155),2));
        time_steps_simulation=t_drying;
        Tin_C=xx.Tin-273.15;
        Tin=xx.Tin;

        for i =2:length(time_steps_HL)
            Tin_C=(1-a)*Tin_C+a*(u.Tinlet_drying-273.15-(a_hl+ug0*...
              ug_exp)*(u.Tinlet_drying-295));              
            Tin(i)=Tin_C+273.15;
        end
        Tin=interp1(time_steps_HL,Tin,time_steps_simulation);          
        
       
        % simulation     
        for i = 1:1:length(time_steps_simulation)-1                    
            % calculation of hT and mT coefficients
            cp_g=p.cp(Tin(i));
            rho_g=101325/8.314/(Tin(i))*0.02897;
            p.h_M=ug0*rho_g/p.zeta/47.2e-6/p.a_V/1e5*d.hM(n_batch);
            p.h_T=cp_g*rho_g*ug0/p.a_V/47.2e-6/p.zeta*d.hT(n_batch);
            
            % the drying model used in this simulator was originally
            % developed for cakes with multiple volatile and non-volatile
            % components of the liquid phase.
            % In this simulator, the number of volatile components is fixed
            % to 1, and the volumetric fraction of non-volatile components
            % of the liquid is null in every node of the grid           
            [~,output]=ode15s(@drying_model,[time_steps_simulation(i) time_steps_simulation(i+1)],x0,options,...
                    xx.number_nodes_drying,p.number_volatile_components,epsL_non_vol,...
                    Tin(i),p.cp_gas_components',p.rho_sol,xx.E,ug0,Pprofile,p.MW_air,...
                    xx.step_grid_drying,p.wl_crit,p.wl_eq,p.coeff_antoine',...
                    p.MW_components,p.h_M,xx.a_V,p.rho_liq_components,p.latent_heat,...
                    p.cp_liq_components,p.cp_s,p.h_T);

            % Extract states from solver output                       
            x0=output(end,:);
            Tgas(i+1,:)=output(end,1:xx.number_nodes_drying);
            Tsolid(i+1,:)=output(end,1+xx.number_nodes_drying:2*xx.number_nodes_drying);
            mass_frG(i+1,:)=output(end,xx.number_nodes_drying*2+1:end-...
                    (p.number_volatile_components)*xx.number_nodes_drying);
            vol_cont_impurities(i+1,:)=reshape(output(end,end-...
                    (p.number_volatile_components)*xx.number_nodes_drying+1:end),...
                    xx.number_nodes_drying,p.number_volatile_components)';                                              
        end
        S_final=vol_cont_impurities/xx.E;   
        
        % impurity content calculation  
        rho_cake=vol_cont_impurities.*rho_liq+(1-xx.E)*p.rho_sol;
        composition=mean(vol_cont_impurities./rho_cake.*p.rho_liq_components,2); 
        
        %% Updates states vector (x) and measurement vector (y)
        % states      
        t_drying=t_drying+t+residual_deliq; % switch to process time
        x.(['pos' num2str(pos)]).Tg=Tgas(end,:);
        x.(['pos' num2str(pos)]).Ts=Tsolid(end,:);
        x.(['pos' num2str(pos)]).Tin=Tin(end);
        x.(['pos' num2str(pos)]).S=S_final(end,:);
        x.(['pos' num2str(pos)]).nodes=xx.nodes_deliq;

        x.(['pos' num2str(pos)]).gas_mass_fr_vect=mass_frG(end,:);

        % measurements
        y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch) ]).t=...
            [y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).t,...
            t_drying(2:end)];
        y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).Tg_bot=...
            [y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).Tg_bot, Tgas(2:end,end)'];
        y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).Ts_bot=...
            [y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).Ts_bot, Tsolid(2:end,end)'];
        y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).w_EtOH_gas=...
            [y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).w_EtOH_gas,...
            mass_frG(2:end,[1:p.number_volatile_components]*xx.number_nodes_drying)'];
        y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).w_EtOH_cake=...
            [y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).w_EtOH_cake,...
            composition(2:end)'];
        y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).S=...
            [y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).S,...
            mean(S_final(2:end,:),2)'];
        y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).Tg_top=...
            [y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).Tg_top, Tin(2:end)];
        y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).Vdryer=...
            [y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).Vdryer, (1000*60*ug0*p.A*ones(1,length(Tin)-1))];
        y.sequence.(['batch_' num2str(n_batch)]).drying_duration=...
            y.sequence.(['batch_' num2str(n_batch)]).drying_duration+Dt;
        
    end
end
    
 