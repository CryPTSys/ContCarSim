function [x,y]=model_filtration(current_time,Dt,p,u,x,y,n_cycle,pos) 

    %% Calculations filtration (cake formation)
    if pos > 1
        x.(['pos' num2str(pos)]).biasV_filt=-y.(['pos' num2str(pos-1)]).(['cycle_' num2str(n_cycle)]).V_filt(end);
        x.(['pos' num2str(pos)]).biast_filt=-y.(['pos' num2str(pos-1)]).(['cycle_' num2str(n_cycle)]).t_filt(end);
    else
        x.(['pos' num2str(pos)]).biasV_filt=0; 
        x.(['pos' num2str(pos)]).biast_filt=0;
    end
    
    xx=x.(['pos' num2str(pos)]);
    initial_liq_mass_fr_vect=xx.liq_mass_fr_vect;  
    t=xx.filtration_time;
    
    % Determine the filtration states as function of time                                       
    t_filt=unique([(p.V_filt_step-round(rem(t,p.V_filt_step),6)+t):p.V_filt_step:Dt+t Dt+t]);%linspace(0,Dt,100);
    
    a_filt= xx.visc_liq.*xx.alpha.*xx.c;  % darcy coeff 1
    b_filt= 2*p.A*xx.visc_liq.*p.Rm;      % darcy coeff 2
    c_filt = -(xx.alpha*xx.visc_liq*xx.c*xx.V_filt^2+2*p.A*xx.visc_liq.*p.Rm*xx.V_filt+...
        2*p.A^2*u.dP*(t_filt-t));% darcy coeff 3
    V_filt= (-b_filt+sqrt(b_filt.^2-4.*a_filt.*c_filt))./(2.*a_filt);     % V(t) - filtrate volume [m^3]
    
    %% Collect outputs in the object filtration_output
    % states
    x.(['pos' num2str(pos)]).filtration_time=x.(['pos' num2str(pos)]).filtration_time+Dt;
   % x.(['pos' num2str(pos)]).t_filt=t_filt;
    x.(['pos' num2str(pos)]).V_filt=V_filt(end);
    x.(['pos' num2str(pos)]).S=1;
%     x.(['pos' num2str(pos)]).dP_media_vacuum=u.dP./(x.pos1.alpha*x.pos1.c*x.pos1.V_filt_final/p.A+p.Rm)*p.Rm;
    
    % measurements
    y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).t_filt=...
        [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).t_filt,...
        t_filt+x.(['pos' num2str(pos)]).biast_filt];
    
    y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).V_filt=...
        [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).V_filt,...
        V_filt+x.(['pos' num2str(pos)]).biasV_filt];  
    
    y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).w_filt=...
        [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).w_filt,...
        repmat(initial_liq_mass_fr_vect,[1,...
        length(t_filt(1:end))])];

end