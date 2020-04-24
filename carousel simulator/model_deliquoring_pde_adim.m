function deliq_output = model_deliquoring_pde_adim(input,duration_deliq,p)
    % Required inputs:
    % input.final_liq_mass_fr_vect = vector liquid phase components mass fractions
    %    [number of components x number of nodes] or [number of components x 1] if uniform
    % input.S_final = vector of initial liquid phase saturation [1 x nodes previous model]
    % p = parameters object
    % Optional inputs:
    % input.nodes = vector of nodes previous model [1 x nodes previous model]
    %   necessary if input.final_liq_mass_fr_vect is not uniform
    
    %% Deliquoring model calculations - solution obtained solving the PDEs model
    t_deliq=0:p.time_step_deliq:duration_deliq;
    
    % Grid discretization
    p.number_nodes=p.number_nodes_deliq;% local variable - round(p.L_cake/p.grid_deliq);
    p.nodes_list=1:p.number_nodes;
    p.step_grid_deliq=p.L_cake/p.number_nodes;
    p.nodes_deliq=linspace(p.step_grid_deliq/2,p.L_cake-p.step_grid_deliq/2,p.number_nodes);
    
    % Initial condition
    initial_liq_mass_fr_vect=input.final_liq_mass_fr_vect;
    if size(initial_liq_mass_fr_vect,2)>1 && (initial_liq_mass_fr_vect(1,1)-initial_liq_mass_fr_vect(1,end))/initial_liq_mass_fr_vect(1)>0.01      
        for i = 1 : size(initial_liq_mass_fr_vect,1)
            initial_liq_mass_fr_vect_rescaled(i,:) = interp1(input.nodes,initial_liq_mass_fr_vect(i,:),p.nodes_deliq);
        end
    else
        initial_liq_mass_fr_vect_rescaled = repmat(initial_liq_mass_fr_vect(:,1),[1,p.number_nodes]);
    end
    S0=input.S_final;
    S0 = interp1(input.nodes,S0,p.nodes_deliq);
     
    % Calculations
    visc_liq=p.visc_liquid_phase_from_mass_fr(mean(initial_liq_mass_fr_vect_rescaled,2));
    dPgdz=p.dP/p.L_cake;
    p.Pgin=101325;
    p.Pgout=101325-p.dP;   
    p.Pg=linspace(p.Pgin-dPgdz*p.step_grid_deliq/2,...
        p.Pgout+dPgdz*p.step_grid_deliq/2,p.number_nodes);
    SR0=(S0-p.S_inf)/(1-p.S_inf); 
    p.scaling=1/(p.step_grid_deliq^3); % an additional scaling factor for pressures
    theta_step=p.k*p.Pb*p.scaling*p.time_step_deliq/(visc_liq*p.L_cake^2*p.E*(1-p.S_inf));
    theta_fin=p.k*p.Pb*p.scaling*duration_deliq/(visc_liq*p.L_cake^2*p.E*(1-p.S_inf));
    
    % Sparsity matrix - consistent with the FVM upwind differencing scheme
%     SparsMatrix=tril(ones(p.number_nodes),1); % lower diagonal 
    options = [];%odeset('JPattern',SparsMatrix);
    
    [~,S_t]=ode45(@deliquoring_pde,0:theta_step:theta_fin,SR0,options,p);
    
    S_t=p.S_inf+S_t*(1-p.S_inf);
    S=trapz(p.nodes_deliq,S_t'/(p.L_cake-p.step_grid_deliq)); % average saturation
    vol_fr_liq_phase=S*p.E;%V_liq_pores_deliq/Volume_cake; % Volumetric fraction of solvent content during deliq
    eq_vol_fr_liq_phase=p.S_inf*p.E;%V_liq_pores_eq/Volume_cake; % Volumetric fraction of solvent content at mech eq
    
    %% Collect outputs in the object deliq_output
    deliq_output.t_deliq=t_deliq;
    deliq_output.eq_vol_fr_liq_phase=eq_vol_fr_liq_phase;
    deliq_output.vol_fr_liq_phase=vol_fr_liq_phase;
    deliq_output.S_avg=S;
    deliq_output.S_final=S_t(end,:);
    deliq_output.nodes=p.nodes_deliq;
    deliq_output.final_liq_mass_fr_vect=initial_liq_mass_fr_vect_rescaled;
end

% HRFVM
function dSdtheta=deliquoring_pde(~,x,p)
    p.step_grid_deliq=p.step_grid_deliq/p.L_cake;
    
    SR=x(1:length(p.nodes_deliq));
    Plin=(p.Pgin-p.Pb.*SR(1).^(-1/p.lambda))/p.Pb/p.scaling;
    Plout=(min(p.Pgout-p.Pb*SR(end)^(-1/p.lambda),p.Pgout))/p.Pb/p.scaling;
    ulin=0; 
    
    % Neglect transient
    Pl=(p.Pg'-p.Pb.*SR.^(-1/p.lambda))/p.Pb/p.scaling;
    
    % Approximate transient before air breakthrough
%     if SR(end)<1-1e-6
%         Pl=(p.Pg'-p.Pb.*SR.^(-1/p.lambda))/p.Pb/p.scaling;
%     else
%         Pl=linspace(p.Pg(1)-p.Pb.*SR(1).^(-1/p.lambda),p.Pg(end),p.number_nodes)/p.Pb/p.scaling;
%         Pl=Pl';
%     end

    Pl=Pl(:);
    kl=SR.^((2+3*p.lambda)/p.lambda);
    
    % HRFVM
    flux_Pl_in=Pl(1);
    fluxes_Pl=fluxes_HRFVM_vL(Pl,p);
    dPldz(1)=(fluxes_Pl(1)-flux_Pl_in)*2/p.step_grid_deliq;
    dPldz(2:p.number_nodes)=(fluxes_Pl(2:p.number_nodes)-fluxes_Pl(1:p.number_nodes-1))/p.step_grid_deliq;
    dPldz(p.number_nodes)=(fluxes_Pl(p.number_nodes)-fluxes_Pl(p.number_nodes-1))*2/p.step_grid_deliq;
    
    % FVM - doesn't work, to be fixed
%     flux_Pl_in=Plin;
%     fluxes_Pl=fluxes_FVM_upwind(Pl,p);
%     fluxes_Pl(1)=0.5*(Pl(1)+Plin);
%     dPldz(1)=(fluxes_Pl(1)-flux_Pl_in)*2/p.step_grid_deliq;
%     dPldz(2:p.number_nodes)=(fluxes_Pl(2:p.number_nodes)-fluxes_Pl(1:p.number_nodes-1))/p.step_grid_deliq;
%     dPldz(p.number_nodes)=(fluxes_Pl(p.number_nodes)-fluxes_Pl(p.number_nodes-1))*2/p.step_grid_deliq;

    ul=-kl.*dPldz';
    
    % HRFVM
    flux_ul_in=ulin;
    fluxes_ul=fluxes_HRFVM_vL(ul,p);
    duldz(1)=(fluxes_ul(1)-flux_ul_in)/p.step_grid_deliq;
    duldz(2:p.number_nodes)=(fluxes_ul(2:p.number_nodes)-fluxes_ul(1:p.number_nodes-1))/p.step_grid_deliq;
    duldz(p.number_nodes)=(fluxes_ul(p.number_nodes)-fluxes_ul(p.number_nodes-1))*2/p.step_grid_deliq;

    % FVM - doesn't work, to be fixed
%     flux_ul_in=ulin;
%     fluxes_ul=fluxes_FVM_upwind(ul,p);
%     fluxes_ul(1)=ul(1);
%     duldz(1)=(fluxes_ul(1)-flux_ul_in)*2/p.step_grid_deliq;
%     duldz(2:p.number_nodes)=(fluxes_ul(2:p.number_nodes)-fluxes_ul(1:p.number_nodes-1))/p.step_grid_deliq;
% %     duldz(p.number_nodes)=(fluxes_ul(p.number_nodes)-fluxes_ul(p.number_nodes-1))*2/p.step_grid_deliq;

    dSdtheta=-duldz';
end

function fluxes=fluxes_HRFVM_vL(x,p)
    r(2:p.number_nodes-1)=(x(2:p.number_nodes-1)-x(1:p.number_nodes-2)+1e-8)./...
        (x(3:p.number_nodes)-x(2:p.number_nodes-1)+1e-8);
    phi_r=((r+abs(r))./(1+abs(r)))'; % the first element is a zero
%     phi_r=[0;(max(max(0,min(1,2*r)),min(r,2)))'];
    fluxes(1)=0.5*(x(1)+x(2));
    fluxes(2:p.number_nodes-1)=x(2:p.number_nodes-1)+0.5*phi_r(2:p.number_nodes-1).*(x(3:p.number_nodes)-x(2:p.number_nodes-1));
    fluxes(p.number_nodes)=x(end);
end

function fluxes=fluxes_FVM_upwind(x,p)
    fluxes(1)=0.5*(x(1)+x(2));
    fluxes(2:p.number_nodes)=x(1:p.number_nodes-1);
end
