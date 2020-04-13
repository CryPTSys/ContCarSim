function washing_output=model_washing_SaturatedCakeApprox(filt_output,deliq_output_pde,p) 
    % This model uses the analytical solution for washing saturated cakes, and then adjusts
    % the resulting solvent content with a suitable correlation Wcorr=f(W,Sin)
    
    S=deliq_output_pde.S_final;
    c0=1; % if you have more than one component, adjust with the densities here and in the output section
    k0=0.15;
    gamma=5;
    H_cake=filt_output.H_cake(end);
    u=p.dP/(p.visc*(p.alpha*p.rho_sol*H_cake*(1-p.E)+p.Rm));
    v=u/p.E;
    Q=u*p.A;
    Wf=Q*p.t_rot/(p.E*p.A*H_cake);
    ReSc=v*p.m1/p.m0/p.Di;
    Dl=p.Di*(sqrt(2)^-1+55.5*ReSc^0.96);
    z_vect=deliq_output_pde.nodes_deliq;
    
    % Washing of saturated cake
    W=Wf;
    m=0;
    for z = z_vect
        m=m+1;
        g=@(w) exp(v*z./(2*Dl)-gamma*(1-w.^(-2))*W-W*v*z./(4*Dl*w.^2)-w.^2*v.*z./(4*W*Dl));
        % with back-flux (longer computation)
        c_sat(m)=c0*(1-0.5*(erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*erfc((1+W)./...
            (2*sqrt(W)).*sqrt(v*z/Dl)))+k0*sqrt(v*z./(pi*W*Dl))*integral(g,1,inf));
        % without back-flux
%         c_sat(m)=c0*(1-0.5*(erfc((1-W)./(2*sqrt(W)).*sqrt(v*z/Dl))+exp(v*z/Dl).*erfc((1+W)./...
%                     (2*sqrt(W)).*sqrt(v*z/Dl))));
    end
% end   
    % Compensation for partially deliquored cake
    Wcorr=Wf+15.1*(1-S).*exp(-1.56*c_sat/c0)-7.4*(1-S.^2).*exp(-1.72*c_sat/c0);
    m=0;
    for z = z_vect
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
    washing_output.xv_mother_liquor=c*p.E;
    washing_output.xv_wash_solvent=(1-c)*p.E;
    washing_output.washing_ratio=Wf;
    washing_output.washing_volume=Q*p.t_rot;
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
