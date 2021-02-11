function [x,y] = model_deliquoring_grad(cycle_time,Dt,p,u,x,y,n_cycle,pos)
    % differential deliquoring model, withOUT species balance: to be used
    % when initial liquid composition is UNIFORM

    %%  Inputs list
    %   cycle_time = cycle timer
    %   Dt = duration of deliquoring step
    %   p = properties object
    %   u = inputs vector
    %   x = states (+additional properties) object
    %   y = measurements vector
    %   n_cycle = cycle number
    %   pos = port in which filtration is occurring   

    %% update equilibrium saturation to current pressure drop
    rho_liq=p.rho_liquid_phase_from_mass_fr(mean(x.(['pos' num2str(pos)]).liq_mass_fr_vect,2));
    x.(['pos' num2str(pos)]).S_inf=trapz(x.(['pos' num2str(pos)]).x,0.155*(1+...
        0.031*p.N_cap_CSD(x.(['pos' num2str(pos)]).x,...
        rho_liq,x.(['pos' num2str(pos)]).E,...
        u.dP,x.(['pos' num2str(pos)]).L_cake).^(-0.49)).*...
        x.(['pos' num2str(pos)]).CSD/x.(['pos' num2str(pos)]).m0);       

    %% if cake is small, use design charts for simulating deliquoring
    if x.(['pos' num2str(pos)]).L_cake <p.min_length_discr
        [x,y] = model_deliquoring(cycle_time,Dt,p,u,x,y,n_cycle,pos);
    else
    %% Deliquoring model calculations - solution obtained solving the PDEs model
        % Initial conditions (creation of local objects with shorter names for conciseness)    
        t=cycle_time;
        xx=x.(['pos' num2str(pos)]);         
        initial_liq_mass_fr_vect=xx.liq_mass_fr_vect;       
        initial_liq_mass_fr_vect_rescaled = ones(size(initial_liq_mass_fr_vect,1),xx.number_nodes_deliq).*initial_liq_mass_fr_vect(:,1);
        S0=xx.S;
        S0 = ones(1,xx.number_nodes_deliq).*S0(:,1);

        % Calculations
        visc_liq=p.visc_liquid_phase_from_mass_fr(298,mean(initial_liq_mass_fr_vect_rescaled,2)); % modify room temperature if needed
        dPgdz=(u.dP-xx.dP_media_vacuum)/xx.L_cake;
        Pgin=101325;
        Pgout=101325-(u.dP-xx.dP_media_vacuum);   
        Pg=linspace(Pgin-dPgdz*xx.step_grid_deliq/2,...
            Pgout+dPgdz*xx.step_grid_deliq/2,xx.number_nodes_deliq);
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
        
        % Saturation and filtrate volume calculation
        S_t=xx.S_inf+S_t*(1-xx.S_inf); % saturation dynamic profile
        V_filt=sum((S0-S_t)*xx.E*xx.step_grid_deliq*p.A,2); % filtrate as function of time [m3]

        %% Updates states vector (x) and measurement vector (y)
        % states
        t_deliq=t+t_deliq;
        V_filt=V_filt';
        x.(['pos' num2str(pos)]).S=S_t(end,:);
        x.(['pos' num2str(pos)]).nodes=xx.nodes_deliq;
        x.(['pos' num2str(pos)]).liq_mass_fr_vect=initial_liq_mass_fr_vect_rescaled;
    
        % measurements
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).t_filt=...
            [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).t_filt...
            t_deliq(2:end)];
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).V_filt=...
            [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).V_filt, ...
            y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).V_filt(end)+...
            V_filt(2:end)];
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).w_filt=...
            [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).w_filt,...
            repmat(initial_liq_mass_fr_vect(:,end),[1,...
            length(t_deliq(2:end))])];
        
    end
end
