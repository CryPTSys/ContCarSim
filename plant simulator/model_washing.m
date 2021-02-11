function [x,y] = model_washing(cycle_time,Dt,p,u,x,y,n_cycle,pos)
    %%  Inputs list
    %   cycle_time - unused
    %   Dt = duration of washing step
    %   p = properties object
    %   x = states (+additional properties) object
    %   u = inputs vector
    %   y = measurements vector
    %   n_cycle = cycle number
    %   pos = port in which filtration is occurring 
    
    %% Washing model calculation - analytical solution
    % Initial conditions (creation of local objects with shorter names for conciseness)    
    xx=x.(['pos' num2str(pos)]);
    S=xx.S; % m3_liq/m3_pores
    initial_liq_mass_fr_vect=xx.liq_mass_fr_vect;
    t=cycle_time;
    
    % Rescale initial saturation and composition profiles to match current grid
    if length(S)>1 
        ws = warning('off','all');
        poly_fit_S=polyfit(xx.nodes,S,4);
        warning(ws)
        S=polyval(poly_fit_S,xx.nodes_washing);
    else
        S=S*ones(1,xx.number_nodes_washing);
    end
    initial_liq_mass_fr_vect_rescaled=zeros(p.number_components,xx.number_nodes_washing);
    if (initial_liq_mass_fr_vect(1,1)-initial_liq_mass_fr_vect(1,end))/initial_liq_mass_fr_vect(1)>0.01
        for i = 1 : size(initial_liq_mass_fr_vect,1)
            initial_liq_mass_fr_vect_rescaled(i,:) = interp1(xx.nodes,initial_liq_mass_fr_vect(i,:),xx.nodes_washing,'pchip');
        end
    else
        initial_liq_mass_fr_vect_rescaled=ones(size(initial_liq_mass_fr_vect,1),...
            xx.number_nodes_washing).*initial_liq_mass_fr_vect(:,1);
    end
    c0=initial_liq_mass_fr_vect_rescaled*p.rho_liquid_phase_from_mass_fr(mean(initial_liq_mass_fr_vect_rescaled,2)); % initial concentration of the species in the liquid phase [kg_i/m3_liq]

    % Washing parameters and pre-calculations
    v=xx.Q_wash/p.A/xx.E;    % wash solvent pore velocity
    t_wash=unique([0 (p.filtration_sampling_time/10-... % wash solvent time vector - starts from 0
        round(rem(t,p.filtration_sampling_time/10),6)):p.filtration_sampling_time/10:Dt Dt]);
    W_vector=(xx.V_wash+xx.Q_wash*t_wash)/(xx.E*p.A*xx.L_cake); % vector of washing ratios corresponding to t_wash
    initial_empty_volume=sum((1-S)*xx.E*xx.step_grid_washing*p.A,2); % volume of empty cake at washing onset        
    ReSc=v*xx.m1/xx.m0./p.Di_liq; % dispersion equation for every species in the liquid phase
    Dl=p.Di_liq.*(sqrt(2)^-1+55.5*ReSc.^0.96); % axial diffusivity coefficients
    p.lambda_ads=(1+p.k_ads*(1-xx.E)/xx.E)^-1;
    
    % Calculate time profiles at outlet for measurements
    liq_mass_fr_vect=[];
    for n = 1:length(W_vector)
        W=W_vector(n);
        liq_mass_fr_vect(:,n)=p.mass_fr_from_conc(species_balance_end_point(c0(:,end),p.c_inlet,Dl,S(end),W,v,p,xx.L_cake,xx));
    end    
    
    % calculate final composition
    c=species_balance(c0,p.c_inlet,Dl,S,W,v,p,xx);
    final_liq_mass_fr_vect=p.mass_fr_from_conc(c); % mass fractions in the liquid phase

    %% Updates states vector (x) and measurement vector (y)
    t_filt=t+t_wash; % time vector of current simulation step - for measurements, starts from cycle time before washing
    V_filt=max(xx.Q_wash*t_wash-initial_empty_volume,0); % volume of wash solvent filtered at current simulation step
    
    % states
    x.(['pos' num2str(pos)]).V_wash=x.(['pos' num2str(pos)]).V_wash+xx.Q_wash*Dt; % volume of filtered wash solvent
    x.(['pos' num2str(pos)]).S=1;
    x.(['pos' num2str(pos)]).nodes=xx.nodes_washing;
    x.(['pos' num2str(pos)]).liq_mass_fr_vect=final_liq_mass_fr_vect;
    x.(['pos' num2str(pos)]).washing_time=x.(['pos' num2str(pos)]).washing_time+Dt; % washing timer
    
    % measurements
    y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).t_filt=...
            [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).t_filt,...
            t_filt(2:end)];
    y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).V_filt=...
            [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).V_filt ...
            V_filt(2:end)+y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).V_filt(end)];
    y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).w_filt=...
            [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).w_filt,...
            liq_mass_fr_vect(:,2:end)];
    
    % if washing hasn't finished, ripristinate composition at washing onset 
    % in the state vector: washing solution works starting from conditions at
    % washing beginning
    if x.(['pos' num2str(pos)]).V_wash<x.(['pos' num2str(pos)]).V_wash_final-1e-10
        x.(['pos' num2str(pos)]).liq_mass_fr_vect=initial_liq_mass_fr_vect_rescaled;
    end
        
end

function c=species_balance(c0,c_inlet,Dl,S,W,v,p,xx)
    
    z = ones(p.number_components,length(xx.nodes_washing)).*xx.nodes_washing;
    Dl= ones(size(Dl,1),size(z,2)).*Dl;
    W_lambda=W*p.lambda_ads;
    c_saturated_cake=c0+(c_inlet-c0).*(0.5*(erfc((z/xx.L_cake-W_lambda)./...
                        (2*sqrt(W_lambda)).*sqrt(v*z./Dl))+exp(v*z./Dl).*...
                        erfc((z/xx.L_cake+W_lambda)./(2*sqrt(W_lambda)).*sqrt(v*z./Dl))));
                    
    % Compensatation for pre-deliquoring - local Wcorr
    Wcorr=W+15.1*(1-S).*exp(-1.56*(c_saturated_cake-c_inlet)./(c0-c_inlet))-...
        7.4*(1-S.^2).*exp(-1.72*(c_saturated_cake-c_inlet)./(c0-c_inlet));
    Wcorr=Wcorr*p.lambda_ads;
    c=c0+(c_inlet-c0).*(0.5*(erfc((z/xx.L_cake-Wcorr)./(2*sqrt(Wcorr)).*...
        sqrt(v*z./Dl))+exp(v*z./Dl).*erfc((z/xx.L_cake+Wcorr)./...
                    (2*sqrt(Wcorr)).*sqrt(v*z./Dl))));
               
end

function c=species_balance_end_point(c0,c_inlet,Dl,S,W,v,p,L_cake,xx)
    
    z = L_cake;
    Dl= ones(size(Dl,1),size(z,2)).*Dl;
    W_lambda=W*p.lambda_ads;
    c_saturated_cake=c0+(c_inlet-c0).*0.5.*(erfc((1-W_lambda)./...
                        (2*sqrt(W_lambda)).*sqrt(v*z./Dl))+exp(v*z./Dl).*...
                        erfc((1+W_lambda)./(2*sqrt(W_lambda)).*sqrt(v*z./Dl)));
%     c=c_saturated_cake;
    % Compensatation for pre-deliquoring - local Wcorr
    S=max(S);
    Wcorr=W+15.1*(1-S).*exp(-1.56*(c_saturated_cake-c_inlet)./(c0-c_inlet))-...
        7.4*(1-S.^2).*exp(-1.72*(c_saturated_cake-c_inlet)./(c0-c_inlet));
    Wcorr=Wcorr*p.lambda_ads;
    c=c0+(c_inlet-c0).*0.5.*(erfc((1-Wcorr)./(2*sqrt(Wcorr)).*...
        sqrt(v*z./Dl))+exp(v*z./Dl).*erfc((1+Wcorr)./...
                    (2*sqrt(Wcorr)).*sqrt(v*z./Dl)));
               
end

