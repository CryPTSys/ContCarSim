function deliq_output = rto_model_deliquoring(input,duration_deliq,p)
    % uses design charts for attaining deliquoring output
    
    % Inputs:
    % input.S_final                     = avg initial liquid phase saturation 
    % p                                 = parameters object

    % Outputs:
    % deliq_output.S_final = final avg cake saturation
    % deliq_output.ThetaP = deliquoring duration in dimensionless time units 

    %% Deliquoring model calculation - solution obtained through design charts
    
    visc_liq=p.visc_liq_components(p.T_room);
    t_deliq=duration_deliq;
    
    % if the cake is not saturated at t=0, calculate corresponding ThetaP0
    SR0=(mean(input.S_final)-p.S_inf)/(1-p.S_inf);
    if  SR0 > .95
        ThetaP0=0;
    elseif SR0 >.3
        ThetaP0=((1-SR0)/1.08/SR0)^(1/.88);
    elseif SR0==0
        ThetaP0=1e5;
    else
        ThetaP0=((1-SR0)/1.46/SR0)^(1/.48);
    end
    
    ThetaP=(t_deliq.*p.k.*(p.dP))./(p.E.*visc_liq*(p.L_cake.^2).*(1-p.S_inf))+ThetaP0; % Dimensionless deliquoring time [-]  

    n_switch=min(max(ThetaP-1.915,0)*100,1);   
    SR=1./(1+1.46.*ThetaP.^0.48).*n_switch+1./(1+1.08.*ThetaP.^0.88).*(1-n_switch);
    S=p.S_inf+SR.*(1-p.S_inf); 

    deliq_output.S_final=S(end);
    deliq_output.ThetaP=ThetaP;
    
end
