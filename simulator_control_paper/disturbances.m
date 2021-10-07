function [cryst_output,d,p]=disturbances(t,cryst_output,p,d,u,n_cycle,flag)
    % Function called to simulate disturbances
    
    %% mesh fouling and cleaning-in-place routine   
    p.Rm=d.resistances(n_cycle,:)';

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
