function [x,y]=model_filtration(~,Dt,p,u,x,y,n_batch,pos) 
    %%  Inputs list
    %   cycle_time - unused
    %   Dt = duration of filtration step
    %   p = properties object
    %   x = states (+additional properties) object
    %   y = measurements vector
    %   n_cycle = cycle number
    %   pos = port in which filtration is occurring 
         
    %% Simulation   
    % Initial conditions (creation of local objects with shorter names for conciseness)
    xx=x.(['pos' num2str(pos)]); % copy of state (+ additional properties) object
    t=xx.filtration_time;
    
    % Create filtration time vector accordingly to sampling time, but without exceeding filtration duration                                    
    t_filt=unique([(p.filtration_sampling_time-round(rem(t,p.filtration_sampling_time),6)+t):p.filtration_sampling_time:Dt+t Dt+t]);%linspace(0,Dt,100);
    
    a_filt= xx.visc_liq.*xx.alpha.*xx.c;    % darcy coeff 1
    b_filt= 2*p.A*xx.visc_liq.*p.Rm(pos);        % darcy coeff 2
    c_filt = -(xx.alpha*xx.visc_liq*xx.c*xx.V_filt^2+2*p.A*xx.visc_liq.*p.Rm(pos)*xx.V_filt+...
        2*p.A^2*u.dP*(t_filt-t));           % darcy coeff 3
    V_filt= (-b_filt+sqrt(b_filt.^2-4.*a_filt.*c_filt))./(2.*a_filt);   % V(t) - filtrate volume [m^3]
    
    %% Updates states vector (x) and measurement vector (y)
    % States
    x.(['pos' num2str(pos)]).filtration_time=x.(['pos' num2str(pos)]).filtration_time+Dt;
    x.(['pos' num2str(pos)]).V_filt=V_filt(end);
    x.(['pos' num2str(pos)]).S=1;
    
    % If filtration started in a previous station, consider only filtrate
    % volume and filtration time carried out in current station in
    % measurement vector y
    if pos > 1
        x.(['pos' num2str(pos)]).V_filt0=-y.(['pos' num2str(pos-1)]).(['cycle_' num2str(n_batch)]).V_filt(end);
        x.(['pos' num2str(pos)]).t_filt0=-y.(['pos' num2str(pos-1)]).(['cycle_' num2str(n_batch)]).t_filt(end);
    else
        x.(['pos' num2str(pos)]).V_filt0=0; 
        x.(['pos' num2str(pos)]).t_filt0=0;
    end
     
    % Measurements
    y.(['pos' num2str(pos)]).(['cycle_' num2str(n_batch)]).t_filt=...
        [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_batch)]).t_filt,...
        t_filt+x.(['pos' num2str(pos)]).t_filt0];

    y.(['pos' num2str(pos)]).(['cycle_' num2str(n_batch)]).V_filt=...
        [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_batch)]).V_filt,...
        V_filt+x.(['pos' num2str(pos)]).V_filt0];  


    %% Check if filtration finished
    if abs(x.(['pos' num2str(pos)]).filtration_time-x.(['pos' num2str(pos)]).filtration_duration)<1e-2
        x.(['pos' num2str(pos)]).filtration_finished=1;
    end

end

