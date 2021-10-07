function y = continuous_outputs(t,p,x,y,n_rotation)
    % This function is called at the end of each cycle to concatenate the
    % measurements collected during the last cycle with the previous one.
    % The output vectors, contained in y.cont_sign, contain continuous
    % samplings, mimicking a real sensor installed on a physical unit

    %% Vectors creation
    if n_rotation == 1
        y.cont_sign.pos1.t_filt=[0 t];
        y.cont_sign.pos1.V_filt=[0 0];
    end
    if n_rotation == 2
        y.cont_sign.pos2.t_filt=[0 t];
        y.cont_sign.pos2.V_filt=[0 0];
    end
    if n_rotation == 3
        y.cont_sign.pos3.t_filt=[0 t];
        y.cont_sign.pos3.V_filt=[0 0];
    end    
    if n_rotation == 4
        y.cont_sign.pos4.t_filt=[0 t];
        y.cont_sign.pos4.V_filt=[0 0];
        y.cont_sign.pos4.w_filt=zeros(p.number_components,2);
        y.cont_sign.pos4.wg=zeros(p.number_volatile_components,2);
        y.cont_sign.pos4.t_drying=[0 t];
        y.cont_sign.pos4.Tg=[298 298];
        y.cont_sign.pos4.Ts=[298 298];
        y.cont_sign.pos4.Tin=[298 298];
    end
    
    %% Sensors Position 1
    if n_rotation > 0        
        y.cont_sign.pos1.t_filt=[y.cont_sign.pos1.t_filt y.cont_sign.pos1.t_filt(end)+...
            y.pos1.(['cycle_' num2str(n_rotation)]).t_filt];
        y.cont_sign.pos1.V_filt=[y.cont_sign.pos1.V_filt,...
            y.pos1.(['cycle_' num2str(n_rotation)]).V_filt];           
        y.cont_sign.pos1_4.bias=y.cont_sign.pos1_4.V(end);
    end
    %% Sensors Position 2
    if n_rotation > 1
        y.cont_sign.pos2.t_filt=[y.cont_sign.pos2.t_filt y.cont_sign.pos2.t_filt(end)+...
            y.pos2.(['cycle_' num2str(n_rotation-1)]).t_filt];

        y.cont_sign.pos2.V_filt=[y.cont_sign.pos2.V_filt,...
            y.pos2.(['cycle_' num2str(n_rotation-1)]).V_filt];

    end
    %% Sensors Position 3
    if n_rotation > 2
        y.cont_sign.pos3.t_filt=[y.cont_sign.pos3.t_filt y.cont_sign.pos3.t_filt(end)+...
            y.pos3.(['cycle_' num2str(n_rotation-2)]).t_filt];

        y.cont_sign.pos3.V_filt=[y.cont_sign.pos3.V_filt,...
            y.pos3.(['cycle_' num2str(n_rotation-2)]).V_filt];

    end
    
    %% Sensors Position 4
    if n_rotation > 3
        y.cont_sign.pos4.t_filt=[y.cont_sign.pos4.t_filt y.cont_sign.pos4.t_filt(end)+...
            y.pos4.(['cycle_' num2str(n_rotation-3)]).t_filt];

        y.cont_sign.pos4.V_filt=[y.cont_sign.pos4.V_filt,...
            y.pos4.(['cycle_' num2str(n_rotation-3)]).V_filt];
        
        
        y.cont_sign.pos4.t_drying=[y.cont_sign.pos4.t_drying,...
            y.cont_sign.pos4.t_drying(end)+...
            y.pos4.(['cycle_' num2str(n_rotation-3)]).t_drying];
        
        y.cont_sign.pos4.Tg=[y.cont_sign.pos4.Tg,...
            y.pos4.(['cycle_' num2str(n_rotation-3)]).Tg];
        
        y.cont_sign.pos4.Tin=[y.cont_sign.pos4.Tin,...
            y.pos4.(['cycle_' num2str(n_rotation-3)]).Tin];
        
        y.cont_sign.pos4.Ts=[y.cont_sign.pos4.Ts,...
            y.pos4.(['cycle_' num2str(n_rotation-3)]).Ts];
        
        y.cont_sign.pos4.wg=[y.cont_sign.pos4.wg,...
            y.pos4.(['cycle_' num2str(n_rotation-3)]).wg];
        
        %% calculation of composition of discharged cake (mass fractions)       
        eps_l=mean(x.pos4.E.*x.pos4.S);       
        rho_liq=p.rho_liq_components;
        rho_cake=eps_l*rho_liq+(1-x.pos4.E)*p.rho_sol;        
        y.final_composition(n_rotation-3,:)=...
               mean(rho_liq.*x.pos4.E.*...
               x.pos4.S/rho_cake,2);          
    end    
end       