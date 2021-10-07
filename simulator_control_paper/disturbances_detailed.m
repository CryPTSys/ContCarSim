function [cryst_output,d,p]=disturbances(t,cryst_output,p,d,u,n_cycle,flag)
    % Function called to simulate disturbances
    
    %% mesh fouling and cleaning-in-place routine
    load resistances
    
    a=[2e9 2e9 2e9 2e9]';
    p.n_cycle_cip=p.n_cycle_cip+1;
    p.Rm=a+p.Rm;

    p.ports_working=1:4;

    n_fouled=sum(p.Rm>=1.4e10);
    if n_fouled > 0 && p.n_cycle_cip > 4
        p.n_cycle_cip=1;                      
    end
    if p.n_cycle_cip < 4 
       p.Rm(1:p.n_cycle_cip)=p.Rm(1:p.n_cycle_cip)-a(1:p.n_cycle_cip);
       p.ports_working=(p.n_cycle_cip+1):4;
    elseif p.n_cycle_cip == 4
       p.ports_working=1;
       p.Rm=3e9+rand(4,1)*1e9;            
    elseif p.n_cycle_cip < 8
       p.Rm(p.n_cycle_cip-3:end)=p.Rm(p.n_cycle_cip-3:end)-a(p.n_cycle_cip-3:end);
       p.ports_working=1:(p.n_cycle_cip-3);
    end
    
    %% gaussian disturbances
    d.c_slurry_dist=d.c_slurry(n_cycle);
    d.V_slurry_dist=d.V_slurry(n_cycle);     
    d.E_dist=d.E(n_cycle);
    d.alpha_dist=1+(1-d.E_dist)/d.E_dist^3; 
    
    %% slurry concentration ramp
    if flag == 1            
        if t>300
            cryst_output.conc_MSMPR=cryst_output.conc_MSMPR*(1+(t-300)/60*0.02);
        end        
    end
    
end