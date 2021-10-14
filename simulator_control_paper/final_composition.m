function y = final_composition(p,x,y,n_rotation)
    % calculation of composition of discharged cake (mass fractions)      
    
    if n_rotation > 3
        if p.ports_working(4)==1     
            eps_l=mean(x.pos4.E.*x.pos4.S);       
            rho_liq=p.rho_liq_components;
            rho_cake=eps_l*rho_liq+(1-x.pos4.E)*p.rho_sol;        
            y.final_composition(n_rotation-3)=...
                   mean(rho_liq.*x.pos4.E.*...
                   x.pos4.S/rho_cake,2); 
        else
            y.final_composition(n_rotation-3)=0;
        end
    end    
end