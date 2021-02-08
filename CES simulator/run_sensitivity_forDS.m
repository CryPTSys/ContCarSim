% Sensitivity analysis - effect of decision variables on 
rng default

clear,close all %,clc

% CPPs grid definition
t_rot_vector=10:10:220; % [s]
V_slurry_vector=(.5:0.5:10)*1e-6; % [mL]
c_slurry=50;

% initialization
it_number=0;
load CSD
CSD=CSD/100;
CSD=CSD/sum(CSD);
cryst_output.x=x;
cryst_output.CSD=CSD;

tic
% calculation of probability for every grid point
for i =1:length(t_rot_vector)
    for j=1:length(V_slurry_vector)
        it_number=it_number+1;
        
        tot_runs=400;
        
        % vectors pre-allocation
        design_vector(:,it_number)=[t_rot_vector(i),V_slurry_vector(j)];
        impurity_vector=zeros(1,tot_runs);
        hm_dist=zeros(1,tot_runs);
        Rm_dist=zeros(1,tot_runs);
        E_dist=zeros(1,tot_runs);
        alpha_dist=zeros(1,tot_runs);
        V_slurry_dist=zeros(1,tot_runs);
        c_slurry_dist=zeros(1,tot_runs);
        
        % model parameters distribution for uncertainty
        hm_dist= 1+randn(1,tot_runs)*0.05;
        Rm_dist=randi([3E9,3e10],1,tot_runs); 
        E_dist= 1+randn(1,tot_runs)*0.05;
        alpha_dist=1+randn(1,tot_runs)*0.05;
        V_slurry_dist=1+randn(1,tot_runs)*0.03;
        c_slurry_dist=1+randn(1,tot_runs)*0.03;
        
        % calculation of probability at current grid point
        for MC = 1:tot_runs           
            impurity_vector(MC)=run_sensitivity(t_rot_vector(i), V_slurry_vector(j),...
            hm_dist(MC), Rm_dist(MC), E_dist(MC), alpha_dist(MC), V_slurry_dist(MC),...
            c_slurry_dist(MC),cryst_output,c_slurry);
        end
        probability(it_number)=sum(impurity_vector<0.005)/tot_runs;
        save('sens_c_slurry50','design_vector','probability')
    end
    
end
toc

function output=run_sensitivity(t_rot, V_slurry, hm_dist, Rm_dist, E_dist, alpha_dist, ...
    V_slurry_dist, c_slurry_dist,cryst_output,c_slurry)
    %% Setup decision variables
    p=carousel_parameters_class;

    p.t_rot=t_rot; % s
    p.V_slurry=V_slurry;
    p.W=1; % Washing ratio
    p.dP=5e4; % Pressure drop filtration, washing and deliquoring [Pa]
    p.dP_drying=p.dP; % Pressure drop drying [Pa]
    p.Tinlet_drying=70+273.15; % Drying gas temperature [K]
    cryst_output.liq_mass_fr_vect=[0.95 0 0.05]'; % 99% mother liquor, 1% impurity
    cryst_output.T=298;

    % set discretization options
    p.number_nodes_washing=50; % number of nodes washing (analytical solution)
    p.min_length_discr=2e-4;   % grid spacing for deliquoring and drying
    p = carousel_parameters(cryst_output,p);   
    
    % uncertainty
    p.E=p.E*E_dist;
    p.alpha=p.alpha*alpha_dist;
    p.Rm=Rm_dist;
    p.h_M=p.h_M*hm_dist;
    p.V_slurry=V_slurry*V_slurry_dist;
    cryst_output.conc_MSMPR=c_slurry*c_slurry_dist;
    
    %% Simulation  
    output=carousel_simulator_funct_sens(cryst_output,p);
    
end