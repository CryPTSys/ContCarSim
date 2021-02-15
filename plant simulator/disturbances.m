function [cryst_output,p]=disturbances(t,cryst_output,p,u_ss,flag)
        % Function called to simulate disturbances
        % flag = 1: ramp in p.Rm (filter mesh resistance)
        % flag = 2: ramp in impurity content in inlet slurry

        if t > u_ss.t_rot*6 && t < u_ss.t_rot*75
            if flag == 1
                p.Rm= 2e9+(t-6*u_ss.t_rot)*2e9*.001;
            elseif flag == 2
                cryst_output.liq_mass_fr_vect(3)=min(0.01+(t-45)*0.01*0.01,0.05);
                cryst_output.liq_mass_fr_vect(1)=1-cryst_output.liq_mass_fr_vect(3);
            end
        else
            p.Rm=2e9;
        end
        
        % store disturbance profile
        p.Rm_vector=[p.Rm_vector p.Rm];
        p.time_vector=[p.time_vector t];
end