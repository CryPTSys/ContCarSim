% Dryer pressure profile

function f = MomBal_vectorized(P,p,x)
    
    % States matrix
    x_j=zeros(p.n_nodes,5);
    x=[x, P];
    for j = 1:p.n_nodes
        x_j(j,:)=x(j:p.n_nodes:end);
    end
    
    ug_j=p.ug(x_j);
    rho_g_j=p.rho_g(x_j);
    eps_g_j=p.eps_g(x_j);
    
    deltaP=[P(2:end) p.Pg_outlet*2]-[P(1:end-1) P(end)*2];
    deltaP=deltaP';
         
    f = deltaP/p.step_size+150*p.mu.*ug_j./p.d_p.*(1-...
            eps_g_j.^2)./eps_g_j.^3+1.75*rho_g_j.*...
            ug_j.^2/p.d_p.*(1-eps_g_j)./eps_g_j.^3;    

end