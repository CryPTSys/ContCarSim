% simulate carousel operation for a given set of decision variables and inputs

function output=run_simulation()
    %% Setup 
    p=carousel_parameters_class; % create parameters object
    
    % set decision variables and inputs
    p.t_rot=100; % s
    p.V_slurry=3e-6;
    p.W=1; % Washing ratio
    p.wash_solvent_mass_fr=[0 1 0]'; % mass fractions - components 1-3 
    p.dP=5e4; % Pressure drop filtration, washing and deliquoring [Pa]
    p.dP_drying=p.dP; % Pressure drop drying [Pa]
    p.Tinlet_drying=70+273.15; % Drying gas temperature [K]
    cryst_output.conc_MSMPR=50; % kg/m3   
    cryst_output.liq_mass_fr_vect=[0.95 0 0.05]'; % 95% mother liquor, 5% impurity
    cryst_output.T=298;

    % load CSD
    load CSD
    cryst_output.x=x;
    cryst_output.CSD=CSD/100/(pi/6)./x.^4; % number based distribution
    cryst_output.x_perc=x_perc;
    cryst_output.CSD_perc=CSD_perc/100; % volume percentage distribution
 
    % set discretization options
    p.number_nodes_washing=50; % number of nodes washing (analytical solution)
    p.min_length_discr=2e-4;   % grid spacing for deliquoring and drying
    
    % in previous versions, you could specify the time step for
    % obtaining dynamic profiles of the system states. Doesn't work in
    % current version, but can be easily changed (see comments in
    % models functions)
%         p.time_step_deliq=p.t_rot;
%         p.time_step_filt=p.t_rot;
%         p.time_step_drying=p.t_rot;
    
    % calculate remaining parameters
    p = carousel_parameters(cryst_output,p);   
    
    tic
    %% Simulation  
    output=carousel_simulator(cryst_output,p);
    toc
