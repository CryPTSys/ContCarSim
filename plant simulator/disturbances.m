function [cryst_output,p,u]=disturbances(t,cryst_output,p,u,flag)
        
        if t > 120*6 && t < 120*75
            if flag == 1
                p.Rm= 3e9+(t-6*120)*3e9*.005;
            elseif flag == 2
                cryst_output.liq_mass_fr_vect(3)=min(0.01+(t-45)*0.01*0.01,0.05);
                cryst_output.liq_mass_fr_vect(1)=1-cryst_output.liq_mass_fr_vect(3);
            end
        else
            p.Rm=3e9;
        end
        p.Rm_vector=[p.Rm_vector p.Rm];
        p.time_vector=[p.time_vector t];
end