function washing_output=model_washing(deliq_output,p) 
    % This model uses the analytical solution for washing saturated cakes, and then adjusts
    % the resulting solvent content with a suitable correlation Wcorr(z)=f(W,S(z))
    
    p.t_washing=0:0.1:p.t_rot;
    p.number_nodes=p.number_nodes_washing;%round(p.L_cake/p.grid_washing);
    p.nodes_list=1:p.number_nodes;
    p.step_grid_washing=p.L_cake/p.number_nodes;
    p.nodes_washing=linspace(p.step_grid_washing/2,p.L_cake-p.step_grid_washing/2,p.number_nodes);
        
    S=deliq_output.S_final;
    if abs(deliq_output.nodes_deliq(4)-p.nodes_washing(4))>1e-10
        ws = warning('off','all');
        poly_fit_S=polyfit(deliq_output.nodes_deliq,S,4);
        warning(ws)
        S=polyval(poly_fit_S,p.nodes_washing);
    end

    c0=1; % if you have more than one component, adjust with the densities here and in the output section
    c_inlet=0;
    k0=0.15;
    gamma=5;
    u=p.dP/(p.visc*(p.alpha*p.rho_sol*p.L_cake*(1-p.E)+p.Rm));
    v=u/p.E;
    Q=u*p.A;    
    W=p.W;%Q*p.t_rot/(p.E*p.A*p.L_cake);
    ReSc=v*p.m1/p.m0/p.Di;
    Dl=p.Di*(sqrt(2)^-1+55.5*ReSc^0.96);
    
    % Washing of saturated cake
    m=0;
    for z = p.nodes_washing
        m=m+1;
        g=@(w) exp(v*z./(2*Dl)-gamma*(1-w.^(-2))*W-W*v*z./(4*Dl*w.^2)-w.^2*v.*z./(4*W*Dl));
        % with back-flux (longer computation)
        c_sat(m)=c0*(1-0.5*(erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*erfc((1+W)./...
            (2*sqrt(W)).*sqrt(v*z/Dl)))+k0*sqrt(v*z./(pi*W*Dl))*integral(g,1,inf));
        % without back-flux
%             c_sat(n,m)=c_inlet+(c0-c_inlet)*(1-0.5*(erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*erfc((1+W)./...
%                         (2*sqrt(W)).*sqrt(v*z/Dl))));
    end

    % Compensation for partially deliquored cake
    Wcorr=W+15.1*(1-S).*exp(-1.56*c_sat/c0)-7.4*(1-S.^2).*exp(-1.72*c_sat/c0);
    m=0;
    for z = p.nodes_washing
        m=m+1;
        W=Wcorr(m);
        g=@(w) exp(v*z./(2*Dl)-gamma*(1-w.^(-2))*W-W*v*z./(4*Dl*w.^2)-w.^2*v.*z./(4*W*Dl));
        % with back-flux (longer computation)
        c(m)=c0*(1-0.5*(erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*erfc((1+W)./...
            (2*sqrt(W)).*sqrt(v*z/Dl)))+k0*sqrt(v*z./(pi*W*Dl))*integral(g,1,inf));
        % without back-flux
%         c(m)=c0*(1-0.5*(my_erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*my_erfc((1+W)./...
%                     (2*sqrt(W)).*sqrt(v*z/Dl))));
    end
    
    % Compensation for partially deliquored cake - avg pre-deliquoring approximation
    W=mean(Wcorr);
    m=0;
    g=@(w) exp(v*z./(2*Dl)-gamma*(1-w.^(-2))*W-W*v*z./(4*Dl*w.^2)-w.^2*v.*z./(4*W*Dl));
    for z = p.nodes_washing
        m=m+1;
        % with back-flux (longer computation)
        c_uniform(m)=c0*(1-0.5*(erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*erfc((1+W)./...
            (2*sqrt(W)).*sqrt(v*z/Dl)))+k0*sqrt(v*z./(pi*W*Dl))*integral(g,1,inf));
        % without back-flux
%         c(m)=c0*(1-0.5*(my_erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*my_erfc((1+W)./...
%                     (2*sqrt(W)).*sqrt(v*z/Dl))));
    end
    
    % Compensation for partially deliquored cake - max pre-deliquoring approximation
    W=max(Wcorr);
    m=0;
    g=@(w) exp(v*z./(2*Dl)-gamma*(1-w.^(-2))*W-W*v*z./(4*Dl*w.^2)-w.^2*v.*z./(4*W*Dl));
    for z = p.nodes_washing
        m=m+1;
        % with back-flux (longer computation)
        c_max(m)=c0*(1-0.5*(erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*erfc((1+W)./...
            (2*sqrt(W)).*sqrt(v*z/Dl)))+k0*sqrt(v*z./(pi*W*Dl))*integral(g,1,inf));
        % without back-flux
%         c(m)=c0*(1-0.5*(my_erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*my_erfc((1+W)./...
%                     (2*sqrt(W)).*sqrt(v*z/Dl))));
    end
    
    % Compensation for partially deliquored cake - min pre-deliquoring approximation
    W=min(Wcorr);
    m=0;
    g=@(w) exp(v*z./(2*Dl)-gamma*(1-w.^(-2))*W-W*v*z./(4*Dl*w.^2)-w.^2*v.*z./(4*W*Dl));
    for z = p.nodes_washing
        m=m+1;
        % with back-flux (longer computation)
        c_min(m)=c0*(1-0.5*(erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*erfc((1+W)./...
            (2*sqrt(W)).*sqrt(v*z/Dl)))+k0*sqrt(v*z./(pi*W*Dl))*integral(g,1,inf));
        % without back-flux
%         c(m)=c0*(1-0.5*(my_erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*my_erfc((1+W)./...
%                     (2*sqrt(W)).*sqrt(v*z/Dl))));
    end
    
%     washing_output.c=c;
    washing_output.xv_mother_liquor=c*p.E;
    washing_output.xv_wash_solvent=(1-c)*p.E;
    washing_output.washing_volume=Q*p.t_rot;
    plot(p.nodes_washing,c_sat(end,:)),hold on,plot(p.nodes_washing,c(end,:)),plot(p.nodes_washing,c_uniform(end,:)),
    plot(p.nodes_washing,c_max(end,:)),plot(p.nodes_washing,c_min(end,:))
    legend('without pre-deliquoring','W=Wcorr(z)','W=mean(Wcorr(z))','W=max(Wcorr(z))','W=min(Wcorr(z))')
end

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
