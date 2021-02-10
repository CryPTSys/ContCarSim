function [x,y] = model_deliquoring(t,Dt,p,u,x,y,n_cycle,pos)
    % uses design charts if starting from saturated cake and if no
    % extrapolation is needed (10% tolerance), else solves PDE system

    %% Deliquoring model calculation - solution obtained through design charts
    xx=x.(['pos' num2str(pos)]);
        
    % Initial condition 
    initial_liq_mass_fr_vect=xx.liq_mass_fr_vect;       
    S0=xx.S;
  
    % Calculations    
    visc_liq=p.visc_liquid_phase_from_mass_fr(298,mean(initial_liq_mass_fr_vect,2)); % modify room temperature if needed  
    t_deliq=unique([0 (p.V_filt_step-round(rem(t,p.V_filt_step),6)):p.V_filt_step:Dt Dt]);%linspace(0,Dt,100);
    
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
    
    ThetaP=(t_deliq.*xx.k.*(u.dP))./(xx.E.*visc_liq*(xx.L_cake.^2).*(1-xx.S_inf))+ThetaP0; % Dimensionless deliquoring time [-]  
    
%     if ThetaP>204*1.1 || input.S<1*0.99 % solve PDE system    
%         deliq_output=model_deliquoring_pde_adim(input,duration_deliq,p);
%         deliq_output.ThetaP=ThetaP;
%         deliq_output.extrapolation=1;
%     else % use design charts
    n_switch=min(max(ThetaP-1.915,0)*100,1);   
    SR=1./(1+1.46.*ThetaP.^0.48).*n_switch+1./(1+1.08.*ThetaP.^0.88).*(1-n_switch);
    S=xx.S_inf+SR.*(1-xx.S_inf); 
    
    t_filt=t_deliq+t;
    V_filt=(mean(S0)-S)*p.A*xx.L_cake*xx.E;
    x.(['pos' num2str(pos)]).S=S(end);
    x.(['pos' num2str(pos)]).nodes=[];
    x.(['pos' num2str(pos)]).liq_mass_fr_vect=repmat(initial_liq_mass_fr_vect(:,end),[1,length(ThetaP)]);

    % measurements
    y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).t_filt=...
        [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).t_filt ...
        t_filt(2:end)];
    y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).V_filt=...
        [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).V_filt, ...
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).V_filt(end)+...
        V_filt(2:end)];
    y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).w_filt=...
        [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).w_filt,...
        repmat(x.(['pos' num2str(pos)]).liq_mass_fr_vect(:,end),[1 length(ThetaP)-1])];

        

end
