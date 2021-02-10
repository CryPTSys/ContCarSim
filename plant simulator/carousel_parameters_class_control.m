classdef carousel_parameters_class_control
    properties
        names_components
        Rm_vector
        time_vector
        control_time
        error_t_rot
        V_filt_step
        drying_time_step
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
        filtration_sampling_time
        drying_sampling_time
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
        h_T
        epsL_non_vol
        step_grid_deliq
        nodes_deliq
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
        wash_solvent_mass_fr
        k_ads
    end 
    methods 
        function output=visc_liq_components(obj,T)  % Liquid viscosity [Pa s] - T in [K]
            output=10.^(obj.visc_liq_coeff(:,1)+obj.visc_liq_coeff(:,2)/T+...
                obj.visc_liq_coeff(:,3)*T+obj.visc_liq_coeff(:,4)*T^2)*1e-3;  
        end
        function output = mass_fr_from_conc(obj,c) % converts concentrations into mass fractions - liquid phase
            output=c./sum(c); % sum(c) is the density of the liquid phase
        end
        function output = mol_fr_from_mass_fr(obj,w) % converts mass fractions into mole fractions - liquid phase
            output=w./obj.MW_components./sum(w./obj.MW_components);
        end
        function output = rho_liquid_phase_from_mass_fr(obj,w) % calculates density of liquid phase from mass fractions
            output=1./sum(w./obj.rho_liq_components); 
        end
        function output = visc_liquid_phase_from_mass_fr(obj,T,w) % calculates viscosity of liquid phase from mass fractions
            output=exp(sum(obj.mol_fr_from_mass_fr(w).*log(visc_liq_components(obj,T))));
        end
        function output = cp_liquid_phase_from_mass_fr(obj,w) % calculates specific heat of liquid phase from mass fractions
            output=sum(w.*obj.cp_liq_components);
        end
        function output = conc_from_mass_fr(obj,w) % calculates concentrations from mass fractions - liquid phase
            output= w.*obj.rho_liquid_phase_from_mass_fr(w);
        end
              
        function output=alpha_CSD(obj,dp,E) % calculates specific cake resistance of cake with particles of size dp and porosity E - Kozeny–Carman eq
            output=180*(1-E)./(E^3*dp.^2*obj.rho_sol);
        end
        function output=N_cap_CSD(obj,x,rho_liq,E,dP,L_cake) % capillary number - Tarleton and Wakeman, 2007
            output=(E^3*(x).^2*(rho_liq*9.81.*L_cake+dP))/((1-E)^2*L_cake.*obj.surf_t);
        end
        function output=pb_CSD(obj,x,E) % threshold pressure - Tarleton and Wakeman, 2007
            output=4.6*(1-E)*obj.surf_t./(E*x);
        end
    end
end
