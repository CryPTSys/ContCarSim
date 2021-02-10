function y = continuous_outputs(t,p,x,y,n_rotation)
    
    % vectors creation
    if n_rotation == 1
        y.cont_sign.pos1.t_filt=[0 t];
        y.cont_sign.pos1.V_filt=[0 0];
        y.cont_sign.pos1.w_filt=zeros(p.number_components,2);
        
    end
    if n_rotation == 2
        y.cont_sign.pos2.t_filt=[0 t];
        y.cont_sign.pos2.V_filt=[0 0];
        y.cont_sign.pos2.w_filt=zeros(p.number_components,2);
    end
    if n_rotation == 3
        y.cont_sign.pos3.t_filt=[0 t];
        y.cont_sign.pos3.V_filt=[0 0];
        y.cont_sign.pos3.w_filt=zeros(p.number_components,2);
    end
    
    if n_rotation == 4
        y.cont_sign.pos4.t_filt=[0 t];
        y.cont_sign.pos4.V_filt=[0 0];
        y.cont_sign.pos4.w_filt=zeros(p.number_components,2);
        y.cont_sign.pos4.wg=zeros(p.number_volatile_components,2);
        y.cont_sign.pos4.t_drying=[0 t];
        y.cont_sign.pos4.Tg=[298 298];
    end
    
    % adding new measurements
    if n_rotation > 0
        
        y.cont_sign.pos1.t_filt=[y.cont_sign.pos1.t_filt y.cont_sign.pos1.t_filt(end)+...
            y.pos1.(['cycle_' num2str(n_rotation)]).t_filt];
        y.cont_sign.pos1.V_filt=[y.cont_sign.pos1.V_filt,...
            y.pos1.(['cycle_' num2str(n_rotation)]).V_filt];
        
        y.cont_sign.pos1.w_filt=[y.cont_sign.pos1.w_filt,...
            y.pos1.(['cycle_' num2str(n_rotation)]).w_filt];
        
        
           
    end
    if n_rotation > 1
        y.cont_sign.pos2.t_filt=[y.cont_sign.pos2.t_filt y.cont_sign.pos2.t_filt(end)+...
            y.pos2.(['cycle_' num2str(n_rotation-1)]).t_filt];

        y.cont_sign.pos2.V_filt=[y.cont_sign.pos2.V_filt,...
            y.pos2.(['cycle_' num2str(n_rotation-1)]).V_filt];
        
        y.cont_sign.pos2.w_filt=[y.cont_sign.pos2.w_filt,...
            y.pos2.(['cycle_' num2str(n_rotation-1)]).w_filt];
    end
    if n_rotation > 2
        y.cont_sign.pos3.t_filt=[y.cont_sign.pos3.t_filt y.cont_sign.pos3.t_filt(end)+...
            y.pos3.(['cycle_' num2str(n_rotation-2)]).t_filt];

        y.cont_sign.pos3.V_filt=[y.cont_sign.pos3.V_filt,...
            y.pos3.(['cycle_' num2str(n_rotation-2)]).V_filt];
        
        y.cont_sign.pos3.w_filt=[y.cont_sign.pos3.w_filt,...
            y.pos3.(['cycle_' num2str(n_rotation-2)]).w_filt];
    end
    if n_rotation > 3


        y.cont_sign.pos4.t_filt=[y.cont_sign.pos4.t_filt y.cont_sign.pos4.t_filt(end)+...
            y.pos4.(['cycle_' num2str(n_rotation-3)]).t_filt];

        y.cont_sign.pos4.V_filt=[y.cont_sign.pos4.V_filt,...
            y.pos4.(['cycle_' num2str(n_rotation-3)]).V_filt];
        
        y.cont_sign.pos4.w_filt=[y.cont_sign.pos4.w_filt,...
            y.pos4.(['cycle_' num2str(n_rotation-3)]).w_filt];
        
        y.cont_sign.pos4.t_drying=[y.cont_sign.pos4.t_drying,...
            y.cont_sign.pos4.t_drying(end)+...
            y.pos4.(['cycle_' num2str(n_rotation-3)]).t_drying];
        
        y.cont_sign.pos4.Tg=[y.cont_sign.pos4.Tg,...
            y.pos4.(['cycle_' num2str(n_rotation-3)]).Tg];
        
        y.cont_sign.pos4.wg=[y.cont_sign.pos4.wg,...
            y.pos4.(['cycle_' num2str(n_rotation-3)]).wg];
        
        % final composition
        
        eps_l=mean(x.pos4.E.*x.pos4.S);
        
        mass_fr_vector=x.pos4.liq_mass_fr_vect;
        rho_liq=p.rho_liquid_phase_from_mass_fr(mean(mass_fr_vector,2));
        rho_cake=eps_l*rho_liq+(1-x.pos4.E)*p.rho_sol;
        rhoL=p.rho_liquid_phase_from_mass_fr(mass_fr_vector);
        y.final_composition(n_rotation-3,:)=...
               mean(mass_fr_vector.*rhoL.*x.pos4.E.*...
               x.pos4.S/rho_cake,2);
           

    end    
end       