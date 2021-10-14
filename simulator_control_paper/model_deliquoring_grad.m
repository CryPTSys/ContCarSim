function [x,y] = model_deliquoring_grad(batch_time,Dt,p,u,x,y,n_batch,pos)
    % differential deliquoring model, withOUT species balance: to be used
    % when initial liquid composition is UNIFORM

    %%  Inputs list
    %   batch_time = cycle timer
    %   Dt = duration of deliquoring step
    %   p = properties object
    %   u = inputs vector
    %   x = states (+additional properties) object
    %   y = measurements vector
    %   n_cycle = cycle number
    %   pos = port in which filtration is occurring   

    %% update equilibrium saturation to current pressure drop
    rho_liq=p.rho_liq_components;
    x.(['pos' num2str(pos)]).S_inf=0.085;

    %% if cake is small, use design charts for simulating deliquoring
    if x.(['pos' num2str(pos)]).L_cake <p.min_length_discr
        [x,y] = model_deliquoring(batch_time,Dt,p,u,x,y,n_batch,pos);
    else
    %% Deliquoring model calculations - solution obtained solving the PDEs model
        % Initial conditions (creation of local objects with shorter names for conciseness)    
        t=batch_time;
        xx=x.(['pos' num2str(pos)]);                    
        S0=xx.S;
        
        if size(S0,2)>1 && abs(S0(1,1)-S0(1,end))/S0(1)>1e-5   
            S0 = interp1(xx.nodes,S0,xx.nodes_deliq);
        else
            S0 = ones(1,xx.number_nodes_deliq).*S0(:,1);
        end
        
        % Calculations
        visc_liq=p.visc_liq_components(298); % modify room temperature if needed
        dPgdz=(u.dP-xx.dP_mesh)/xx.L_cake;
        % pressure filtration
        Pgin=101325+u.dP-dPgdz*xx.step_grid_deliq/2; % pressure first node
        Pgout=101325-dPgdz*xx.step_grid_deliq/2;            % pressure last node
        Pg=linspace(Pgin,Pgout,xx.number_nodes_deliq);
        SR0=(S0-xx.S_inf)/(1-xx.S_inf); 
        scaling=1/(xx.step_grid_deliq^4); % an additional scaling factor for pressures
        
        % time vector
        t_deliq=unique([0 (p.filtration_sampling_time-round(rem(t,p.filtration_sampling_time),6)):p.filtration_sampling_time:Dt Dt]);
        theta_interval=t_deliq*xx.k*xx.Pb*scaling/(visc_liq*xx.L_cake^2*xx.E*(1-xx.S_inf)); % dimensionless time
        
        % Sparsity matrix - consistent with the FVM upwind differencing scheme
        SparsMatrix=triu(tril(ones(xx.number_nodes_deliq),1),-1); % lower diagonal + upper diagonal
        options = odeset('JPattern',SparsMatrix);
        
        xx.step_grid_deliq=xx.step_grid_deliq/xx.L_cake; % step grid becomes dimensionless
       
        % solver call
        [~,S_t]=ode15s(@deliquoring_grad_mex,theta_interval,SR0,options,xx.number_nodes_deliq,...
            xx.step_grid_deliq, xx.Pb, p.lambda, scaling, Pgin, Pgout, Pg);
        
        xx.step_grid_deliq=xx.step_grid_deliq*xx.L_cake; % step grid has dimension again
        
        % if t_deliq has only two elements, ode15s considers more time steps
        % the next if loop is necessary for having S_t consistent with t_deliq
        if length(t_deliq)==2
            S_t=[S_t(1,:); S_t(end,:)];
        end
        
        % Saturation, filtrate volume and residual impurities calculation
        S_t=xx.S_inf+S_t*(1-xx.S_inf); % saturation dynamic profile
        V_filt=sum((S0-S_t)*xx.E*xx.step_grid_deliq*p.A,2); % filtrate as function of time [m3]
        vol_cont_impurities=S_t*xx.E;
        rho_cake=vol_cont_impurities.*rho_liq+(1-xx.E)*p.rho_sol;
        composition=mean(vol_cont_impurities./rho_cake.*p.rho_liq_components,2); 
        
        %% Updates states vector (x) and measurement vector (y)
        % states
        t_deliq=t+t_deliq;
        V_filt=V_filt';
        x.(['pos' num2str(pos)]).S=S_t(end,:);
        x.(['pos' num2str(pos)]).nodes=xx.nodes_deliq;
    
        % measurements
        y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).t=...
            [y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).t...
            t_deliq(2:end)];
        y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).m_filt=...
            [y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).m_filt, ...
            y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).m_filt(end)+...
            V_filt(2:end)*p.rho_liq_components];
        y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).w_EtOH_cake=...
            [y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).w_EtOH_cake,...
            composition(2:end)'];
        y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).S=...
            [y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).S,...
            mean(S_t(2:end,:),2)'];
        y.sequence.(['batch_' num2str(n_batch)]).deliquoring_duration=...
            y.sequence.(['batch_' num2str(n_batch)]).deliquoring_duration+Dt;
        
    end
end
