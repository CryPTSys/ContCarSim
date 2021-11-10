function [cryst_output,d,p]=disturbances(t,cryst_output,cryst_output_nominal,p,d,u,n_cycle,flag)
    % Function called to simulate disturbances
    
    %% mesh fouling and cleaning-in-place routine   
    p.Rm=d.resistances(n_cycle,:)'; 
    
    %% slurry concentration ramp change
    if flag == 1            
        if t>300 && t < 1500
            d.c_slurry(n_cycle)=d.c_slurry(n_cycle)*(1+(t-300)/60*0.02);
        elseif t >= 1500
            d.c_slurry(n_cycle)=d.c_slurry(n_cycle)*(1+(1500-300)/60*0.02);
        end
    end
    
    %% cake resistance step change
    if flag == 2
        if t>300
            d.alpha(n_cycle)=d.alpha(n_cycle)*2;     
        end
    end
    
    %% gaussian disturbances
    d.c_slurry_dist=d.c_slurry(n_cycle);
    d.V_slurry_dist=d.V_slurry(n_cycle);     
    d.E_dist=d.E(n_cycle);
    d.alpha_dist=d.alpha(n_cycle);
    
    %% Store nominal slurry concentration 
    cryst_output.conc_slurry_vector(end+1)=cryst_output.conc_slurry;
end