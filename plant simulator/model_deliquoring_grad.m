function [x,y] = model_deliquoring_grad(t,Dt,p,u,x,y,n_cycle,pos)
    
    if x.(['pos' num2str(pos)]).L_cake <p.min_length_discr
        [x,y] = model_deliquoring(t,Dt,p,u,x,y,n_cycle,pos);
    else
        
        %% Deliquoring model calculations - solution obtained solving the PDEs model
        xx=x.(['pos' num2str(pos)]);
        
        % Initial condition
        initial_liq_mass_fr_vect=xx.liq_mass_fr_vect;       
        S0=xx.S;
        if size(initial_liq_mass_fr_vect,2)>1 && abs(initial_liq_mass_fr_vect(1,1)-initial_liq_mass_fr_vect(1,end))/initial_liq_mass_fr_vect(1)>1e-5      
%             initial_liq_mass_fr_vect_rescaled=zeros(p.number_components,x.number_nodes);
%             for i = 1 : size(initial_liq_mass_fr_vect,1)
%                 initial_liq_mass_fr_vect_rescaled(i,:) = interp1(input.nodes,initial_liq_mass_fr_vect(i,:),x.nodes_deliq,'pchip');
%             end
        else
            initial_liq_mass_fr_vect_rescaled = ones(size(initial_liq_mass_fr_vect,1),xx.number_nodes_deliq).*initial_liq_mass_fr_vect(:,1);
        end
        if size(S0,2)>1 && abs(S0(1,1)-S0(1,end))/S0(1)>0.01      
%            S0 = interp1(input.nodes,S0,xx.nodes_deliq);
        else
            S0 = ones(1,xx.number_nodes_deliq).*S0(:,1);
        end

        % Calculations
        visc_liq=p.visc_liquid_phase_from_mass_fr(298,mean(initial_liq_mass_fr_vect_rescaled,2)); % modify room temperature if needed
        dPgdz=(u.dP-xx.dP_media_vacuum)/xx.L_cake;
        p.Pgin=101325;
        p.Pgout=101325-(u.dP-xx.dP_media_vacuum);   
        p.Pg=linspace(p.Pgin-dPgdz*xx.step_grid_deliq/2,...
            p.Pgout+dPgdz*xx.step_grid_deliq/2,xx.number_nodes_deliq);
        SR0=(S0-xx.S_inf)/(1-xx.S_inf); 
        p.scaling=1/(xx.step_grid_deliq^4); % an additional scaling factor for pressures
%         theta_step=xx.k*xx.Pb*p.scaling*.3/(visc_liq*xx.L_cake^2*xx.E*(1-xx.S_inf));
%         theta_fin=xx.k*xx.Pb*p.scaling*Dt/(visc_liq*xx.L_cake^2*xx.E*(1-xx.S_inf));
        %t_deliq=unique([0 (p.V_filt_step-rem(t,p.V_filt_step)):p.V_filt_step:Dt Dt]);%linspace(0,Dt,100);
        theta_interval=unique([0 (p.V_filt_step-round(rem(t,p.V_filt_step),6)):p.V_filt_step:Dt Dt])*...
            xx.k*xx.Pb*p.scaling/(visc_liq*xx.L_cake^2*xx.E*(1-xx.S_inf));%linspace(0,Dt,100);
        % Sparsity matrix - consistent with the FVM upwind differencing scheme
    %     SparsMatrix=tril(ones(p.number_nodes),1); % lower diagonal 
        SparsMatrix=triu(tril(ones(xx.number_nodes_deliq),1),-1); % lower diagonal + upper diagonal
        options = odeset('JPattern',SparsMatrix);
        xx.step_grid_deliq=xx.step_grid_deliq/xx.L_cake;
       
        [t_del,S_t]=my_ode15s(@deliquoring_grad_mex,theta_interval,SR0,options,xx.number_nodes_deliq,...
            xx.step_grid_deliq, xx.Pb, p.lambda, p.scaling, p.Pgin, p.Pgout, p.Pg);
        xx.step_grid_deliq=xx.step_grid_deliq*xx.L_cake;
        if length(theta_interval)==2
            t_del=[t_del(1) t_del(end)];
            S_t=[S_t(1,:); S_t(end,:)];
        end
       % Filtrate volume calculation
%         Pl=(repmat(p.Pg(end-1:end),[size(S_t,1),1])-xx.Pb.*S_t(:,end-1:end).^(-1/p.lambda))/xx.Pb/p.scaling;
%         kl=S_t(:,end).^((2+3*p.lambda)/p.lambda);
%         dPldz=Pl(:,end)-Pl(:,end-1);
%         ul=-kl.*dPldz;
%         V_filt=ul*p.A*xx.E*(xx.k*xx.Pb*p.scaling/(visc_liq*xx.L_cake^2*xx.E*(1-xx.S_inf)));
        S_t=xx.S_inf+S_t*(1-xx.S_inf);
        V_filt=sum((S0-S_t)*xx.E*xx.step_grid_deliq*p.A,2);
       
%        S=trapz(p.nodes_deliq,S_t'/(p.L_cake-p.step_grid_deliq)); % average saturation
%        vol_fr_liq_phase=S*xx.E;%V_liq_pores_deliq/Volume_cake; % Volumetric fraction of solvent content during deliq
%        eq_vol_fr_liq_phase=p.S_inf*x.E;%V_liq_pores_eq/Volume_cake; % Volumetric fraction of solvent content at mech eq

        %% Collect outputs in the object deliq_output
        % states
        t_filt=t+t_del'/(xx.k*xx.Pb*p.scaling/(visc_liq*xx.L_cake^2*xx.E*(1-xx.S_inf)));
        V_filt=V_filt';
        x.(['pos' num2str(pos)]).S=S_t(end,:);
        x.(['pos' num2str(pos)]).nodes=xx.nodes_deliq;
        x.(['pos' num2str(pos)]).liq_mass_fr_vect=initial_liq_mass_fr_vect_rescaled;
    
        % measurements
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).t_filt=...
            [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).t_filt...
            t_filt(2:end)];
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).V_filt=...
            [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).V_filt, ...
            y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).V_filt(end)+...
            V_filt(2:end)];
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).w_filt=...
            [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).w_filt,...
            repmat(initial_liq_mass_fr_vect(:,end),[1,...
            length(t_filt(2:end))])];
        
    end
end
