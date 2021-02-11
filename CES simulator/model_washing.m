function washing_output=model_washing(input,p) 
    % This model uses the analytical solution for washing saturated cakes, and then adjusts
    % the resulting solvent content with the correlation Wcorr(z)=f(W,S(z))
    % The active code calculates only the axial profiles at the end of the
    % washing step. Previous versions calculated dynamic profiles
    % (check inactive code in final section, to be debugged, for dynamic profiles)
    %
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
    
    
    %% Initialization    
    if p.W==0 % no washing if washiong ratio is 0
        washing_output=input;
        washing_output.washing_ratio=p.W;
        washing_output.washing_duration=0;
    else   
        % Grid discretization
        number_nodes=p.number_nodes_washing;%round(p.L_cake/p.grid_washing);
        p.step_grid_washing=p.L_cake/number_nodes;
        p.nodes_washing=linspace(p.step_grid_washing/2,p.L_cake-p.step_grid_washing/2,number_nodes);       

        % Get initial conditions from input object
        S=input.S_final; % m3_liq/m3_pores
        initial_liq_mass_fr_vect=input.final_liq_mass_fr_vect;
        
        % Rescale initial saturation on washing grid if previous model grid is different  
        % we use a polynomial interpolation to remove discontinuities in S
        if length(S)>1 
            ws = warning('off','all');
            poly_fit_S=polyfit(input.nodes,S,4);
            warning(ws)
            S=polyval(poly_fit_S,p.nodes_washing);
        else
            S=S*ones(1,number_nodes);
        end
        % Rescale initial composition - if not uniform - on washing grid
        initial_liq_mass_fr_vect_rescaled=zeros(p.number_components,number_nodes);
        if (initial_liq_mass_fr_vect(1,1)-initial_liq_mass_fr_vect(1,end))/initial_liq_mass_fr_vect(1)>0.01
            for i = 1 : size(initial_liq_mass_fr_vect,1)
                initial_liq_mass_fr_vect_rescaled(i,:) = interp1(input.nodes,initial_liq_mass_fr_vect(i,:),p.nodes_washing,'pchip');
            end
        else
            initial_liq_mass_fr_vect_rescaled=ones(size(initial_liq_mass_fr_vect,1),...
                number_nodes).*initial_liq_mass_fr_vect(:,1);
        end
        c0=initial_liq_mass_fr_vect_rescaled*p.rho_liquid_phase_from_mass_fr(mean(initial_liq_mass_fr_vect_rescaled,2)); % initial concentration of the species in the liquid phase [kg_i/m3_liq]

        % Washing parameters and pre-calculations
        visc_liq=p.visc_liquid_phase_from_mass_fr(298,[0 1 0]');  % liquid phase viscosity = viscosity washing liquid
        u=p.dP/(visc_liq*(p.alpha*p.rho_sol*p.L_cake*(1-p.E)+p.Rm)); % superficial velocity
        v=u/p.E; % pore velocity (at the beginning the actual pore velocity is higher, because S<1)
        Q=u*p.A; % Washing liquid flowrate  
        W=p.W;%Q*p.t_rot/(p.E*p.A*p.L_cake);
%         washing_duration=W*p.E*p.A*p.L_cake/Q; 
                                                   
        %% Dispersion equation for every species in the liquid phase
        ReSc=v*p.m1/p.m0./p.Di_liq;
        Dl=p.Di_liq.*(sqrt(2)^-1+55.5*ReSc.^0.96); % axial diffusivity coefficients  
        p.lambda_ads=(1+p.k_ads*(1-p.E)/p.E)^-1;

        % Analytical solution neglecting back-flux
        c=species_balance(c0,p.c_inlet,Dl,S,W,v,p);

        final_liq_mass_fr_vect=p.mass_fr_from_conc(c); % mass fractions in the liquid phase
        
        washing_volume=p.W*p.V_liquid_pores;      
        washing_duration=washing_volume/Q; % mean(S) is added because at the beginning
                                                       % the cake is partially empty  
        %% Collect outputs in the object washing_output
        washing_output.final_liq_mass_fr_vect=final_liq_mass_fr_vect;
        washing_output.washing_volume=washing_volume;
        washing_output.nodes=p.nodes_washing;
        washing_output.washing_ratio=W;
        washing_output.washing_duration=washing_duration;
        washing_output.S_final=ones(1,number_nodes);

    end
end

function c=species_balance(c0,c_inlet,Dl,S,W,v,p)

    z = ones(p.number_components,length(p.nodes_washing)).*p.nodes_washing;
    Dl= ones(size(Dl,1),size(z,2)).*Dl;
    
    W_lambda=W*p.lambda_ads;
    c_saturated_cake=c0+(c_inlet-c0).*0.5.*(erfc((z/p.L_cake-W_lambda)./...
                        (2*sqrt(W_lambda)).*sqrt(v*p.L_cake./Dl))+exp(v*z./Dl).*...
                        erfc((z/p.L_cake+W_lambda)./(2*sqrt(W_lambda)).*sqrt(v*p.L_cake./Dl)));

    % Compensatation for pre-deliquoring - local Wcorr
    S=max(S);
    Wcorr=W+15.1*(1-S).*exp(-1.56*(c_saturated_cake-c_inlet)./(c0-c_inlet))-...
        7.4*(1-S.^2).*exp(-1.72*(c_saturated_cake-c_inlet)./(c0-c_inlet));
    Wcorr=Wcorr*p.lambda_ads;
    c=c0+(c_inlet-c0).*0.5.*((erfc((z/p.L_cake-Wcorr)./(2*sqrt(Wcorr)).*...
        sqrt(v*p.L_cake./Dl))+exp(v*z./Dl).*erfc((z/p.L_cake+Wcorr)./...
                    (2*sqrt(Wcorr)).*sqrt(v*p.L_cake./Dl))));
               
end


%% for getting time profiles, debug the following code
% n=0;
%     for W = linspace(0,100,100)
% %
%         n=n+1;
%         m=0;
%             for z = z_vect
%                 m=m+1;
%                 c0_f(n)=k0*exp(-gamma*W);
%                 g=@(w) exp(v*z./(2*Dl)-gamma*(1-w.^(-2))*W-W*v*z./(4*Dl*w.^2)-w.^2*v.*z./(4*W*Dl));
%                 % without back-flux
%                 % c(n,m)=c_inlet+(c0-c_inlet).*(0.5*(erfc((z/p.L_cake-W)./(2*...
%                 %     sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*erfc((z/p.L_cake+W)./...
% %                     (2*sqrt(W)).*sqrt(v*z/Dl))));
%             
%             end
%     end
%     
