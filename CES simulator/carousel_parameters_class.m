classdef carousel_parameters_class
    properties
        wash_solvent_mass_fr
        Di_gas
        V_slurry
        V_filt_step
        drying_time_step
        names_components
        number_components
        number_volatile_components
        visc_liq_coeff 
        MW_components
        rho_liq_components
        c_inlet
        Di_liq
        cp_liq_components
        vl_crit
        h_M
        vl_eq
        latent_heat
        coeff_antoine
        station_diameter
        Rm
        rho_sol
        visc_gas_phase
        surf_t
        lambda
        k0
        gamma
        A
        E
        m0
        m1
        m2
        m3
        alpha
        k
        k_ads
        a_V
        k_air
        MW_air
        cp_air
        cp_s
        number_nodes_deliq
        number_nodes_filt
        number_nodes_washing
        number_nodes_drying
        time_step_deliq
        time_step_filt
        time_step_drying
        min_length_discr
        t_rot
        W
        dP
        dP_drying
        Tinlet_drying
        L_cake
        Pb
        V_liquid_pores
        m_dry_cake
        dP_media_vacuum
        S_inf
        number_nodes
        nodes_list
        step_grid_drying
        nodes_drying
        dPgdz
        ug
        Pprofile
        cp_l
        rho_l
        repLatHeat
        repRhoLComp
        epsL_non_vol
        step_grid_deliq
        nodes_deliq
        h_T
        Pgin
        Pgout
        Pg
        scaling
        step_grid_washing
        nodes_washing
        lambda_ads
        visc_liq
        x
        CSD
        x_perc
        CSD_perc
    end 
    methods 
        function output=visc_liq_components(obj,T) 
            output=10.^(obj.visc_liq_coeff(:,1)+obj.visc_liq_coeff(:,2)/T+...
                obj.visc_liq_coeff(:,3)*T+obj.visc_liq_coeff(:,4)*T^2)*1e-3;  % Liquid viscosity [Pa s] - T in [K]
        end
        function output = mass_fr_from_conc(obj,c)
            output=c./sum(c); % sum(c) is the density of the liquid phase
        end
        function output = mol_fr_from_mass_fr(obj,w)
            output=w./obj.MW_components./sum(w./obj.MW_components);
        end
        function output = rho_liquid_phase_from_mass_fr(obj,w)
            output=1./sum(w./obj.rho_liq_components); 
        end
        function output = visc_liquid_phase_from_mass_fr(obj,T,w)
            output=exp(sum(obj.mol_fr_from_mass_fr(w).*log(visc_liq_components(obj,T))));
        end
        function output = cp_liquid_phase_from_mass_fr(obj,w)
            output=sum(w.*obj.cp_liq_components);
        end
        function output = conc_from_mass_fr(obj,w)
            output= w.*obj.rho_liquid_phase_from_mass_fr(w);
        end
        
        function output=alpha_CSD(obj,dp)
            output=180*(1-obj.E)./(obj.E^3*dp.^2*obj.rho_sol);
        end
        function output=N_cap_CSD(obj,x,rho_liq)
            output=(obj.E^3*(x).^2*(rho_liq*9.81.*obj.L_cake+obj.dP))/((1-obj.E)^2*obj.L_cake.*obj.surf_t);
        end
        function output=pb_CSD(obj,x)
            output=4.6*(1-obj.E)*obj.surf_t./(obj.E*x);      
        end
    end
end
