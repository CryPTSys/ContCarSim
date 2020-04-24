function washing_output=model_washing(input,p) 
    % This model uses the analytical solution for washing saturated cakes, and then adjusts
    % the resulting solvent content with the correlation Wcorr(z)=f(W,S(z))
    % The active code calculates only the axial profiles at the end of the
    % washing step. See commented code below the main function for
    % computing time-varying profiles
    %
    % Required inputs:
    % input.final_liq_mass_fr_vect = vector liquid phase components mass fractions
    %    [number of components x number of nodes] or [number of components x 1] if uniform
    % input.S_final = vector of initial liquid phase saturation [1 x nodes previous model]
    % p = parameters object
    % Optional inputs:
    % input.nodes = vector of nodes previous model [1 x nodes previous model]
    %   necessary if input.final_liq_mass_fr_vect is not uniform
    
    
    %% Initialization    
    % Grid discretization
    number_nodes=p.number_nodes_washing;%round(p.L_cake/p.grid_washing);
    p.step_grid_washing=p.L_cake/number_nodes;
    p.nodes_washing=linspace(p.step_grid_washing/2,p.L_cake-p.step_grid_washing/2,number_nodes);       
    
    % Get initial conditions from input object
    S=input.S_final; % m3_liq/m3_pores
    initial_liq_mass_fr_vect=input.final_liq_mass_fr_vect;
    % Rescale initial saturation on washing grid if previous model grid is different  
    % we use a polynomial interpolation to remove discontinuities in S
    if abs(number_nodes-length(S))>0    
        ws = warning('off','all');
        poly_fit_S=polyfit(input.nodes,S,4);
        warning(ws)
        S=polyval(poly_fit_S,p.nodes_washing);
    end
    % Rescale initial composition - if not uniform - on washing grid
    if (initial_liq_mass_fr_vect(1,1)-initial_liq_mass_fr_vect(1,end))/initial_liq_mass_fr_vect(1)>0.01
        for i = 1 : size(initial_liq_mass_fr_vect,1)
            initial_liq_mass_fr_vect_rescaled(i,:) = interp1(input.nodes,initial_liq_mass_fr_vect(i,:),p.nodes_washing);
        end
    else
        initial_liq_mass_fr_vect_rescaled=repmat(initial_liq_mass_fr_vect(:,1),[1,number_nodes]);
    end
    c0=initial_liq_mass_fr_vect_rescaled*p.rho_liquid_phase_from_mass_fr(mean(initial_liq_mass_fr_vect_rescaled,2)); % initial concentration of the species in the liquid phase [kg_i/m3_liq]
    
    % Washing parameters and pre-calculations
    visc_liq=p.visc_liq_components(2); % liquid phase viscosity = viscosity washing liquid
    u=p.dP/(visc_liq*(p.alpha*p.rho_sol*p.L_cake*(1-p.E)+p.Rm));
    v=u/p.E;
    Q=u*p.A; % Washing liquid flowrate  
    W=p.W;%Q*p.t_rot/(p.E*p.A*p.L_cake);
    washing_duration=W*p.E*p.A*p.L_cake/Q;

    %% Dispersion equation for every species in the liquid phase
    ReSc=v*p.m1/p.m0./p.Di_liq;
    Dl=p.Di_liq.*(sqrt(2)^-1+55.5*ReSc.^0.96); % axial diffusivity coefficients    
%     for i = 1:length(c0)
%         c(i,:)=species_balance(c0(i),p.c_inlet(i),Dl(i),S,W,v,p);
%     end
    for i = 1:length(c0)
        c=species_balance(c0,p.c_inlet',Dl,S,W,v,p);
    end
    final_liq_mass_fr_vect=p.mass_fr_from_conc(c); % mass fractions in the liquid phase
    
    %% Collect outputs in the object washing_output
    washing_output.vol_fr_mother_liquor=c(1,:)/p.rho_liq_components(1)*p.E; % vol. fraction over the cake
    washing_output.vol_fr_wash_solvent=c(2,:)/p.rho_liq_components(2)*p.E; % vol. fraction over the cake
    washing_output.vol_fr_impurity=c(3,:)/p.rho_liq_components(3)*p.E;
    washing_output.final_liq_mass_fr_vect=final_liq_mass_fr_vect;
    washing_output.washing_volume=p.W*p.V_liquid_pores;
    washing_output.nodes=p.nodes_washing;
    washing_output.washing_ratio=W;
    washing_output.washing_duration=washing_duration;
    washing_output.S_final=ones(1,number_nodes);
end

function c=species_balance(c0,c_inlet,Dl,S,W,v,p)
    
    z = repmat(p.nodes_washing,[length(p.names_components),1]);
    Dl= repmat(Dl',[1,size(z,2)]);
    c_saturated_cake=c_inlet+(c0-c_inlet).*(1-0.5*(erfc((1-W)./...
                        (2*sqrt(W)).*sqrt(v*z./Dl))+exp(v*z./Dl).*erfc((1+W)./...
                        (2*sqrt(W)).*sqrt(v*z./Dl))));
    % Simulate washing neglecting pre-deliquoring
%     m=0;
%     c_saturated_cake=zeros(1,number_nodes);
%     for z = p.nodes_washing
%         m=m+1;
% %         g=@(w) exp(v*z./(2*Dl)-p.gamma*(1-w.^(-2))*W-W*v*z./(4*Dl*w.^2)-w.^2*v.*z./(4*W*Dl));
%         % consider diffusive back-flux (longer computation)
% %         c_saturated_cake(m)=c_inlet+(c0-c_inlet)*(1-0.5*(erfc((1-W)./...
% %             (2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*erfc((1+W)./...
% %             (2*sqrt(W)).*sqrt(v*z/Dl)))+p.k0*sqrt(v*z./(pi*W*Dl))*integral(g,1,inf));
%         % neglect diffusive back-flux
%     end

    % Compensatation for pre-deliquoring - local Wcorr
    Wcorr=W+15.1*(1-S).*exp(-1.56*(c_saturated_cake-c_inlet)./(c0-c_inlet))-...
        7.4*(1-S.^2).*exp(-1.72*(c_saturated_cake-c_inlet)./(c0-c_inlet));
%     c=zeros(1,number_nodes);
%     m=0; 
%     for z = p.nodes_washing
%         m=m+1;
%         W=Wcorr(m);
% %         g=@(w) exp(v*z./(2*Dl)-p.gamma*(1-w.^(-2))*W-W*v*z./(4*Dl*w.^2)-w.^2*v.*z./(4*W*Dl));
%         % consider diffusive back-flux (longer computation)
% %         c(m)=c0*(1-0.5*(erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*erfc((1+W)./...
% %             (2*sqrt(W)).*sqrt(v*z/Dl)))+p.k0*sqrt(v*z./(pi*W*Dl))*integral(g,1,inf));
%     end
    
    % neglect diffusive back-flux
    c=c_inlet+(c0-c_inlet).*(1-0.5*(erfc((1-Wcorr)./(2*sqrt(Wcorr)).*sqrt(v*z./Dl))+exp(v*z./Dl).*erfc((1+Wcorr)./...
                    (2*sqrt(Wcorr)).*sqrt(v*z./Dl))));
end
    %% Uncomment the following part for the comparison with Wcorr calculated 
    %  with different assumptions
    
%     % Compensation for partially deliquored cake - avg pre-deliquoring approximation
%     W=mean(Wcorr);
%     m=0;
%     g=@(w) exp(v*z./(2*Dl)-gamma*(1-w.^(-2))*W-W*v*z./(4*Dl*w.^2)-w.^2*v.*z./(4*W*Dl));
%     for z = p.nodes_washing
%         m=m+1;
%         % with back-flux (longer computation)
%         c_uniform(m)=c0*(1-0.5*(erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*erfc((1+W)./...
%             (2*sqrt(W)).*sqrt(v*z/Dl)))+k0*sqrt(v*z./(pi*W*Dl))*integral(g,1,inf));
%         % without back-flux
% %         c(m)=c0*(1-0.5*(my_erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*my_erfc((1+W)./...
% %                     (2*sqrt(W)).*sqrt(v*z/Dl))));
%     end
%     
%     % Compensation for partially deliquored cake - max pre-deliquoring approximation
%     W=max(Wcorr);
%     m=0;
%     g=@(w) exp(v*z./(2*Dl)-gamma*(1-w.^(-2))*W-W*v*z./(4*Dl*w.^2)-w.^2*v.*z./(4*W*Dl));
%     for z = p.nodes_washing
%         m=m+1;
%         % with back-flux (longer computation)
%         c_max(m)=c0*(1-0.5*(erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*erfc((1+W)./...
%             (2*sqrt(W)).*sqrt(v*z/Dl)))+k0*sqrt(v*z./(pi*W*Dl))*integral(g,1,inf));
%         % without back-flux
% %         c(m)=c0*(1-0.5*(my_erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*my_erfc((1+W)./...
% %                     (2*sqrt(W)).*sqrt(v*z/Dl))));
%     end
%     
%     % Compensation for partially deliquored cake - min pre-deliquoring approximation
%     W=min(Wcorr);
%     m=0;
%     g=@(w) exp(v*z./(2*Dl)-gamma*(1-w.^(-2))*W-W*v*z./(4*Dl*w.^2)-w.^2*v.*z./(4*W*Dl));
%     for z = p.nodes_washing
%         m=m+1;
%         % with back-flux (longer computation)
%         c_min(m)=c0*(1-0.5*(erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*erfc((1+W)./...
%             (2*sqrt(W)).*sqrt(v*z/Dl)))+k0*sqrt(v*z./(pi*W*Dl))*integral(g,1,inf));
%         % without back-flux
% %         c(m)=c0*(1-0.5*(my_erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*my_erfc((1+W)./...
% %                     (2*sqrt(W)).*sqrt(v*z/Dl))));
%     end
%     
%     plot(p.nodes_washing/p.L_cake,c_sat(end,:)),hold on,plot(p.nodes_washing/p.L_cake,c(end,:)),plot(p.nodes_washing/p.L_cake,c_uniform(end,:)),
%     plot(p.nodes_washing/p.L_cake,c_max(end,:)),plot(p.nodes_washing/p.L_cake,c_min(end,:))
%     legend('neglecting pre-deliquoring effect','Wcorr=Wcorr(z)','Wcorr=mean(Wcorr(z))','Wcorr=max(Wcorr(z))','Wcorr=min(Wcorr(z))')
%     xlabel('Cake axial coordinate')
%     ylabel([{'c/c0 @ end of washing [-]'}])
%     set(gca,'fontsize',16,'linewidth',1.3,'xtick',0:0.2:1)    


% function output=my_erfc(x)
%     if abs(erf(x)-1)<1e-3
%         output=1-erf(x);
%     else
%         output=erfc(x);
%     end
% end

%% for getting time profiles, use the following code
% n=0;
%     for W = linspace(0,100,100)
% %
%         n=n+1;
%         m=0;
%             for z = z_vect
%                 m=m+1;
%                 c0_f(n)=k0*exp(-gamma*W);
%                 g=@(w) exp(v*z./(2*Dl)-gamma*(1-w.^(-2))*W-W*v*z./(4*Dl*w.^2)-w.^2*v.*z./(4*W*Dl));
%                 % with back-flux (lengthy computation)
%                 cn(n,m)=c0*(1-0.5*(erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*erfc((1+W)./...
%                     (2*sqrt(W)).*sqrt(v*z/Dl)))+k0*sqrt(v*z./(pi*W*Dl))*integral(g,1,inf));
%                 % without back-flux
%                 % c(n,m)=c0*(1-0.5*(erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*erfc((1+W)./...
% %                     (2*sqrt(W)).*sqrt(v*z/Dl))));
%             
%             end
%     end
%     
