function [cryst_output,d,p]=disturbances(t,cryst_output,cryst_output_nominal,p,d,u,n_cycle,flag)
    % Function called to simulate disturbances
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   disturbances are store in object d. See output section below for
    %   information of the structure.
    %
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
    %   p = parameters object (see parameters.m for field values)
    %           field most relevant for disturbances: p.Rm [4 x 1], where
    %           element i is the resistance of the filter mesh of station i for
    %           the current cycle
    %   d = see output section for structure
    %   u   =   vector of manipulated variables for follwing control interval
    %                       Fields of u:
    %                       - u.t_cycle=cycle duration set-point (s)    MUST BE AN INTEGER
    %                       - u.V_slurry=fed slurry volume set-point (m3)
    %                       - u.P_compr= gauge pressure provided by compressor P101 (Pa)
    %                       - u.Tinlet_drying=drying gas temperature Station 4 set-point (K)
    %   n_cycle   =   cycle counter - number of current cycle 
    %   flag = disturbance scenario flag
    %
    %   Outputs
    %   cryst_output: don't use
    %   d = Disturbance object. Update it with desired value of disturbances for the considered
    %       cycle. Fields:
    %
    %           i) values of parametric disturbances at the considered cycle:
    %                   d.V_slurry_dist = parametric disturbance multiplied by slurry volume set-point
    %                                     to calculate actual slurry volume fed
    %                                     to carousel at considered cycle
    %                   d.c_slurry_dist: = parametric disturbance multiplied by nominal slurry conc. 
    %                                      to calculate actual concentration of slurry fed
    %                                      to carousel at considered cycle
    %                 d.E_dist = = parametric disturbance multiplied by nominal
    %                                porosity to calculate porosity of cake
    %                                formed in Station 1 at considered cycle
    %                 d.alpha_dist = parametric disturbance multiplied by nominal
    %                                porosity to calculate porosity of cake
    %                                formed in Station 1 at considered cycle
    %     
    %           ii)vectors storing default (e.g., normal operating conditions)
    %              parametric disturbances profiles for first 1000 cycles of carousel operation. 
    %              Can be modified to imposed abnormal variations, 
    %              e.g. see lines 85-93 for disturbance scenario 1
    %                   d.c_slurry: [1×1000 double]
    %                   d.V_slurry: [1×1000 double]
    %                   d.E: [1×1000 double]
    %                   d.alpha: [1×1000 double]
    %                   d.hM: [1×1000 double]
    %                   d.hT: [1×1000 double]
    %
    %           iii) fouling and cleaning
    %           d.fouling: ignore, redundant
    %           d.resistances = default mesh resistances for first 1200 cycles [1200×4 double]
    %           d.stations_working = working and non working stations for first 1200 cycles [1200×4 double]
    %                                if elements (i,j) is 1, Station k
    %                                is active during cycle i, otherwise
    %                                it is empty during cycle i.
    % 
    %    p = structure defined in parameters.m
    %                       field most relevant for disturbances: p.Rm [4 x 1], where
    %                       element i is the resistance of the filter mesh of station i for
    %                       the current cycle
    %% mesh fouling and cleaning-in-place routine   
    % always called to assign resistances at current cycle (p.Rm)
    p.Rm=d.resistances(n_cycle,:)'; % [4 x 1], element i is the resistance 
                                    % of the filter mesh of Station i for
                                    % the current cycle. To edit the mesh
                                    % resistances schedule, either modify
                                    % resistances.m (from which
                                    % d.resistances is built), or implement
                                    % fouling laws here
    
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
    
    %% assign parametric disturbances for current cycle from vectors storing parametric disturbance profiles
    d.c_slurry_dist=d.c_slurry(n_cycle);
    d.V_slurry_dist=d.V_slurry(n_cycle);     
    d.E_dist=d.E(n_cycle);
    d.alpha_dist=d.alpha(n_cycle);
    % d.hM and d.hT are directlt called inside model_drying.m
    
    %% Store nominal slurry concentration - don't edit
    cryst_output.conc_slurry_vector(end+1)=cryst_output_nominal.conc_slurry;
end