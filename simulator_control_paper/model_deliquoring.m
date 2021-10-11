function [x,y] = model_deliquoring(t,Dt,p,u,x,y,n_batch,pos)
    %%  Inputs list
    %   batch_time - cycle timer
    %   Dt = duration of deliquoring step
    %   p = properties object
    %   x = states (+additional properties) object
    %   y = measurements vector
    %   n_cycle = cycle number
    %   pos = port in which filtration is occurring  

    
    %% update equilibrium saturation to current pressure drop
    rho_liq=p.rho_liquid_phase_from_mass_fr(mean(x.(['pos' num2str(pos)]).liq_mass_fr_vect,2));
    x.(['pos' num2str(pos)]).S_inf=0.085;% sum(0.155*(1+...
%         0.031*p.N_cap_CSD(x.(['pos' num2str(pos)]).x_perc,rho_liq,...
%         x.(['pos' num2str(pos)]).E,u.dP,x.(['pos' num2str(pos)]).L_cake).^...
%         (-0.49)).*x.(['pos' num2str(pos)]).CSD_perc); 
    
    %% Deliquoring model calculation - solution obtained through design charts
    % Initial conditions (creation of local objects with shorter names for conciseness)    
    xx=x.(['pos' num2str(pos)]);
    initial_liq_mass_fr_vect=xx.liq_mass_fr_vect;       
    S0=xx.S;
  
    % Calculations    
    visc_liq=p.visc_liquid_phase_from_mass_fr(298,mean(initial_liq_mass_fr_vect,2)); % modify room temperature if needed  
    t_deliq=unique([0 (p.filtration_sampling_time-round(rem(t,p.filtration_sampling_time),6)):p.filtration_sampling_time:Dt Dt]);%linspace(0,Dt,100);
    
    % if the cake is not saturated at t=0, calculate corresponding ThetaP0
    SR0=(mean(S0)-xx.S_inf)/(1-xx.S_inf);
    if  SR0 > .95
        ThetaP0=0;
    elseif SR0 >.3
        ThetaP0=((1-SR0)/1.08/SR0)^(1/.88);
    elseif SR0==0
        ThetaP0=1e5;
    else
        ThetaP0=((1-SR0)/1.46/SR0)^(1/.48);
    end
    
    % Dimensionless deliquoring time [-]  
    ThetaP=(t_deliq.*xx.k.*(u.dP))./(xx.E.*visc_liq*(xx.L_cake.^2).*(1-xx.S_inf))+ThetaP0; 
    
    % Design charts equation
    n_switch=min(max(ThetaP-1.915,0)*100,1);   
    SR=1./(1+1.46.*ThetaP.^0.48).*n_switch+1./(1+1.08.*ThetaP.^0.88).*(1-n_switch);
    S=xx.S_inf+SR.*(1-xx.S_inf); 
    
    % impurities
    vol_cont_impurities=S_t*xx.E;
    rho_cake=vol_cont_impurities.*rho_liq+(1-xx.E)*p.rho_sol;
    composition=mean(vol_cont_impurities./rho_cake.*p.rho_liq_components,2); 
    
    %% Updates states vector (x) and measurement vector (y)
    % states
    t_filt=t_deliq+t;
    V_filt=(mean(S0)-S)*p.A*xx.L_cake*xx.E;
    x.(['pos' num2str(pos)]).S=S(end);
    x.(['pos' num2str(pos)]).nodes=[];
    x.(['pos' num2str(pos)]).liq_mass_fr_vect=repmat(initial_liq_mass_fr_vect(:,end),[1,length(ThetaP)]);
    
    % measurements
    y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).t=...
        [y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).t ...
        t_filt(2:end)];
    y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).m_filt=...
        [y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).m_filt, ...
        y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).m_filt(end)+...
        V_filt(2:end)*p.rho_liq_components];
    y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).w_EtOH_cake=...
        [y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).w_EtOH_cake,...
        composition(2:end)'];
    y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).w_EtOH_cake=...
        [y.(['pos' num2str(pos)]).(['batch_' num2str(n_batch)]).w_EtOH_cake,...
        S(2:end)'];
    
end
