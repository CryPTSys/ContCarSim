function deliq_output = model_deliquoring_species_grad(input,duration_deliq,p)

    % differential deliquoring model, with species balance: to be used
    % when initial liquid composition is not necessarily uniform

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


    if p.L_cake <p.min_length_discr 
        deliq_output=model_deliquoring(input,duration_deliq,p);     
        
    else
        %% Deliquoring model calculations - solution obtained solving the PDEs model
        % Grid discretization variables
        p.number_nodes=max(p.number_nodes_deliq,2);% local variable - round(p.L_cake/p.grid_deliq);
        p.nodes_list=1:p.number_nodes;
        p.step_grid_deliq=p.L_cake/p.number_nodes;
        p.nodes_deliq=linspace(p.step_grid_deliq/2,p.L_cake-p.step_grid_deliq/2,p.number_nodes);

        % Initial conditions
        initial_liq_mass_fr_vect=input.final_liq_mass_fr_vect;      
        S0=input.S_final;
        % rescale input saturation profile to current cake discretization grid
        if size(initial_liq_mass_fr_vect,2)>1 && abs(initial_liq_mass_fr_vect(1,1)-initial_liq_mass_fr_vect(1,end))/initial_liq_mass_fr_vect(1)>0.01      
            initial_liq_mass_fr_vect_rescaled=zeros(p.number_components,p.number_nodes);
            for i = 1 : size(initial_liq_mass_fr_vect,1)
                initial_liq_mass_fr_vect_rescaled(i,:) = interp1(input.nodes,initial_liq_mass_fr_vect(i,:),p.nodes_deliq,'pchip');
            end
        else
            initial_liq_mass_fr_vect_rescaled = ones(size(initial_liq_mass_fr_vect,1),...
                p.number_nodes).*initial_liq_mass_fr_vect(:,1);
        end
        if size(S0,2)>1 && abs(S0(1,1)-S0(1,end))/S0(1)>0.01      
            S0 = interp1(input.nodes,S0,p.nodes_deliq);
        else
            S0 = ones(1,p.number_nodes).*S0(:,1);
        end

        % Calculations
        p.visc_liq=p.visc_liquid_phase_from_mass_fr(298,mean(initial_liq_mass_fr_vect_rescaled,2)); % modify room temperature if needed
        dPgdz=(p.dP-p.dP_media_vacuum)/p.L_cake;
        p.Pgin=101325;
        p.Pgout=101325-(p.dP-p.dP_media_vacuum);   
        p.Pg=linspace(p.Pgin-dPgdz*p.step_grid_deliq/2,...
            p.Pgout+dPgdz*p.step_grid_deliq/2,p.number_nodes);
        SR0=(S0-p.S_inf)/(1-p.S_inf); 
        p.scaling=1/(p.step_grid_deliq^4); % an additional scaling factor for pressures
%         theta_step=p.k*p.Pb*p.scaling*p.time_step_deliq/(p.visc_liq*p.L_cake^2*p.E*(1-p.S_inf));
        theta_fin=p.k*p.Pb*p.scaling*duration_deliq/(p.visc_liq*p.L_cake^2*p.E*(1-p.S_inf));

        % Sparsity matrix - consistent with the FVM upwind differencing scheme
        SparsMatrix=triu(tril(ones(p.number_nodes),1),-1); % lower diagonal 
    %     SparsMatrix=triu(tril(ones(p.number_nodes)),-1);
        options = odeset('JPattern',SparsMatrix);
        
        p.step_grid_deliq=p.step_grid_deliq/p.L_cake; %  grid step becomes dimensionless
        
        % solver call
        [t_deliq,x]=ode23(@deliquoring_species_grad_mex,[0 theta_fin],[SR0,...
            reshape(initial_liq_mass_fr_vect_rescaled',[1,...
            p.number_nodes*p.number_components])],...
            options, p.number_nodes, p.step_grid_deliq,...
                      p.Pb, p.lambda, p.scaling, p.Pgin, p.Pgout,...
                      p.Pg, p.number_components, p.S_inf);
 
        % output arrays preparation
        p.step_grid_deliq=p.step_grid_deliq*p.L_cake; % grid step has dimension again

        S_t=x(:,1:p.number_nodes); % saturation time and shape profiles
        final_liq_mass_fr_vect_rescaled= reshape(x(end,p.number_nodes+1:end),...
            p.number_nodes,p.number_components)';

        S_t=p.S_inf+S_t*(1-p.S_inf); % absolute saturation time and space profiles (not used in sensitivity)
        S=trapz(p.nodes_deliq,S_t(end,:)'/(p.L_cake-p.step_grid_deliq)); % average saturation
        vol_fr_liq_phase=S*p.E; % Volumetric fraction of solvent content during deliq
        eq_vol_fr_liq_phase=p.S_inf*p.E; % Volumetric fraction of solvent content at mech eq

        %% Collect outputs in the object deliq_output
        deliq_output.t_deliq=t_deliq;
        deliq_output.S_t=S_t; % absolute saturation time and space profiles (not used in sensitivity)
        deliq_output.eq_vol_fr_liq_phase=eq_vol_fr_liq_phase;
        deliq_output.vol_fr_liq_phase=vol_fr_liq_phase;
        deliq_output.S_avg=S;
        deliq_output.S_final=S_t(end,:);
        deliq_output.nodes=p.nodes_deliq;
        deliq_output.final_liq_mass_fr_vect=final_liq_mass_fr_vect_rescaled;
        
                 
    end
end
