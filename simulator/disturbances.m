function [cryst_output,d,p]=disturbances(t,cryst_output,cryst_output_nominal,p,d,u,n_cycle,flag)
    % Function called to simulate disturbance scenarios
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Variability sources characterizing a disturbance scenario are stored
    %   in object d. See output section below for information on the structure.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Inputs
    %   t      =   timer started at process onset (s)
    %   cryst_output: don't use
    %   cryst_output_nominal  = object containing nominal feed
    %                        conditions. Fields:
    %                        - cryst_output_nominal.conc_slurry= nominal
    %                                 slurry concentration in feed (kg/m3)
    %                        - cryst_output_nominal.x = Crystal size
    %                                 distribution – particles diameters (m)
    %                        - cryst_output_nominal.CSD - Volumetric crystal size distribution – percentage
    %                        - cryst_output_nominal.CSD_perc - Volumetric crystal size distribution – percentage
    %                        - cryst_output_nominal.T - slurry temperature (= room temperature) (K)  
    %   p = parameters object (see parameters.m for complete fields list)
    %           field most relevant for disturbance scenarios: p.Rm [4 x 1], where
    %           element i is the resistance of the filter mesh of station i for
    %           the current cycle (variability sources 1-4)
    %   d = see output section for structure
    %   u   =   vector of set-points of operating variables for following control interval
    %           Fields of u:
    %           - u.t_cycle=cycle duration set-point (s)    MUST BE AN INTEGER
    %           - u.V_slurry=fed slurry volume set-point (m3)
    %           - u.P_compr= gauge pressure provided by compressor P101 (Pa)
    %           - u.Tinlet_drying=drying gas temperature Station 4 set-point (K)
    %   n_cycle   =   cycle counter - number of current cycle 
    %   flag = disturbance scenario flag
    %
    %   Outputs
    %   cryst_output: don't use/modify
    %   d = Disturbance object. Contains current values, profiles and other
    %       information for variability sources 1-10
    %           Fields:
    %
    %           i) fouling and cleaning: related to variability sources 1-4
    %           Modify d.resistances for modifying variability sources 1-4
    %           inter-cycle profiles (edit d.stations_working consistently).
    %               d.fouling: ignore, redundant
    %               d.resistances = inter-cycle mesh resistance profiles loaded from resistances.m
    %                               With default resistances.m, d.resistances contains 
    %                               mesh resistances for first 1200 cycles [1200×4 double]
    %                               At cycle i, variability sources 1-4
    %                               assume the value of d.resistances(i,:)
    %               d.stations_working = working and non working stations for first 1200 cycles [1200×4 double]
    %                                if elements (i,j) is 1, Station k
    %                                is active during cycle i, otherwise
    %                                it is empty during cycle i.
    %                                If the resistances pattern in
    %                                resistance.m is modified,
    %                                d.stations_working should be modified
    %                                coherently.
    %     
    %           ii) vectors storing default (e.g., normal operating conditions)
    %              variability sources 5-10 values for first 1000 cycles of carousel operation. 
    %              Modify for modifying variability sources 5-10
    %              inter-cycle profiles.
    %              e.g. see lines 89-96 for disturbance scenario 1
    %                   d.c_slurry: [1×1000 double]
    %                   d.hM: [1×1000 double]
    %                   d.hT: [1×1000 double]
    %                   d.V_slurry: [1×1000 double]
    %                   d.E: [1×1000 double]
    %                   d.alpha: [1×1000 double]
    %
    %           iii) values of variability sources 5,8,9,10 at the considered cycle
    %               (automatically assigned from the respective profile)
    %                   d.c_slurry_dist: = parameter multiplied by nominal slurry conc. 
    %                                      to calculate actual concentration of slurry fed
    %                                      to carousel at considered cycle
    %                   d.V_slurry_dist = parameter multiplied by slurry volume set-point
    %                                     to calculate actual slurry volume fed
    %                                     to carousel at considered cycle
    %                   d.E_dist = = parameter multiplied by nominal
    %                                porosity to calculate porosity of cake
    %                                formed in Station 1 at considered cycle
    %                   d.alpha_dist = parameter multiplied by nominal
    %                                porosity to calculate porosity of cake
    %                                formed in Station 1 at considered cycle
    %
    %    p = structure defined in parameters.m
    %               field most relevant for disturbance scenarios: p.Rm [4 x 1], where
    %               element i is the resistance of the filter mesh of station i for
    %               the current cycle (variability sources 1-4)    
    %% slurry concentration ramp change
    % disturbance scenario 1
    if flag == 1            
        if t>300 && t < 1500
            d.c_slurry(n_cycle)=d.c_slurry(n_cycle)*(1+(t-300)/60*0.02);
        elseif t >= 1500
            d.c_slurry(n_cycle)=d.c_slurry(n_cycle)*(1+(1500-300)/60*0.02);
        end
    end
    
    %% cake resistance step change
    % disturbance scenario 2
    if flag == 2
        if t>300
            d.alpha(n_cycle)=d.alpha(n_cycle)*2;     
        end
    end
    
    %% assign variability sources 1-4 at current cycle (p.Rm) from d.resistances
    % don't edit: modify d.resistances if you want to modify variability
    %             sources 1-4 profiles
    p.Rm=d.resistances(n_cycle,:)'; % [4 x 1], element i is the resistance 
                                    % of the filter mesh of Station i for
                                    % the current cycle. To edit the mesh
                                    % resistances schedule, either modify
                                    % resistances.m (from which
                                    % d.resistances is built), or implement
                                    % fouling laws here. Please edit
                                    % d.stations_working consistently
    %% assign variability sources 5,8,9,10 at current cycle from vectors storing their inter-cycle profiles
    % don't edit: modify the vectors on the RHS if you want to modify the
    %             a variability source profile
    d.c_slurry_dist=d.c_slurry(n_cycle);
    d.V_slurry_dist=d.V_slurry(n_cycle);     
    d.E_dist=d.E(n_cycle);
    d.alpha_dist=d.alpha(n_cycle);
    % d.hM (variability source 6) and d.hT (variability source 7) are directly 
    % called inside model_drying.m: directly modify d.hM(n_cycle) and/or
    % d.hT(n_cycle) for modifying hM_dist and hT_dist
    
    %% Store nominal slurry concentration 
    % don't edit
    cryst_output.conc_slurry_vector(end+1)=cryst_output_nominal.conc_slurry;
end