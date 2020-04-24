% Carousel simulator
% F. Destro, v1: March27, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
close all

%% Load data from file, in the future these data will be provided by cryst model
load CSD_daniel
CSD=CSD(end,:);
cryst_output.x=x;
cryst_output.CSD=CSD;
cryst_output.conc_MSMPR=11; % kg/m3   
cryst_output.flowrate_MSMPR=1.7e-7; % m3/s
cryst_output.liq_mass_fr_vect=[0.99 0 0.01]'; % 99% mother liquor, 1% impurity
clear x, clear CSD

%% Inputs and decision variables
p.t_rot=250; % s
p.dP=5e4; % Pressure drop filtration, washing and deliquoring [Pa]
p.W=5; % Washing ratio
p.dP_drying=101325; % Pressure drop drying [Pa]
p.Tinlet_drying=70+273.15; % Drying gas temperature [K]
t=0:0.01:p.t_rot;
%% Solution and discretization options
p.number_nodes_deliq=50;
p.number_nodes_filt=p.number_nodes_deliq; % used for constructing output vectors, not for computations
p.number_nodes_washing=100;
p.number_nodes_drying=30;
p.time_step_deliq = .1; % Time step for calculating saturation profile during deliquoring [s]
p.time_step_filt = .1;  % Time step for calculating saturation profile during filtration [s]
p.time_step_drying = .1;
%% Parameters of the model and cake properties calculation
p = carousel_parameters(cryst_output,p);

%% Simulation section

%% Position 2: filtration (+ pre-deliquoring)
% Filtration 
[filt_output,p]=model_filtration(cryst_output,p);

% Pre-deliquoring
residual_time_pos2=p.t_rot-filt_output.filtration_duration; % Duration of pre-deliquoring

if residual_time_pos2 > 0 % pre-deliquoring occurs in position 2
    % Solve with design charts 
    deliq_output_pos2_charts=model_deliquoring_design_charts(residual_time_pos2,p);
    % Solve integrating the PDEs
    output_pos2=model_deliquoring_pde_adim(filt_output,residual_time_pos2,p);   
    residual_filtration_duration=0; % starting time of washing in position 3 
    
elseif residual_time_pos2 < 0 % filtration doesn't finish in position 2 and goes on in position 3
    residual_filtration_duration=abs(residual_time_pos2);
    output_pos2=filt_output;
end

%% Positions 3 and 4: (filtration +) washing + deliquoring
if residual_filtration_duration <= p.t_rot % washing starts in position 3
    washing_output=model_washing(output_pos2,p);
    residual_time_pos3=p.t_rot-residual_filtration_duration-washing_output.washing_duration; 
    if residual_time_pos3 > 0 % the cake undergoes a deliquoring step after washing in position 3
         % Uncomment the following line to see the effect of deliquoring in Position 3 after washing
         %output_pos3=model_deliquoring_pde_adim(washing_output,residual_time_pos3,p); 
        output_pos4=model_deliquoring_pde_adim(washing_output,residual_time_pos3+p.t_rot,p);
    elseif residual_time_pos3 == 0 % the cake does not undergo deliquoring after washing in position 3
        output_pos4=model_deliquoring_pde_adim(washing_output,p.t_rot,p);  
    elseif p.t_rot-abs(residual_time_pos3)>=0 %the washing solvent hold-up finishes to be filtered in position 4
        output_pos4=model_deliquoring_pde_adim(washing_output,p.t_rot-abs(residual_time_pos3),p);
    end
elseif residual_filtration_duration<2*p.t_rot % filtration finishes in position4, no time for washing in position 3 
    output_pos4=model_deliquoring_pde_adim(filt_output,p.t_rot-(residual_filtration_duration-p.t_rot),p);
else % filtration is not complete when entering the drying step!
    output_pos4=output_pos2;
    residual_filtration_duration=residual_filtration_duration-p.t_rot; % this residual filtration is to be carried out in position 5
end

% %% Position 5: drying
% drying_output=model_drying(output_pos4_pde,p);

% 
% figure
% plot(drying_output.Tgas(1:200:end,:)')
% figure
% plot(drying_output.xv_liquid_phase(1,:)),hold on,plot(drying_output.xv_liquid_phase(end,:))
%% Print a table with the performed steps
fprintf('\nPosition 1\n')
fprintf('%2.1f mL of slurry loaded\n',cryst_output.flowrate_MSMPR*1e6*p.t_rot)
fprintf('\nPosition 2\n')
if residual_time_pos2 > 0 % pre-deliquoring occurs in position 2
    fprintf('Filtration: %2.1f seconds\n',filt_output.filtration_duration)
    fprintf('Pre-deliquoring: %2.1f seconds\n',residual_time_pos2) 
elseif residual_time_pos2 < 0 % filtration doesn't finish in position 2 and goes on in position 3
     fprintf('Filtration: %2.1f seconds\n',p.t_rot)
end
if residual_filtration_duration <= p.t_rot % washing starts in position 3
   if residual_time_pos3 > 0 % the cake undergoes a deliquoring step after washing in position 3
        if residual_filtration_duration == 0 
            fprintf('\nPosition 3\n')
            fprintf('Washing: %2.1f seconds\n',washing_output.washing_duration)
            fprintf('Deliquoring: %2.1f seconds\n',p.t_rot-residual_filtration_duration-washing_output.washing_duration)
        else
            fprintf('\nPosition 3\n')
            fprintf('Filtration: %2.1f seconds\n',residual_filtration_duration)
            fprintf('Washing: %2.1f seconds\n',washing_output.washing_duration)
            fprintf('Deliquoring: %2.1f seconds\n',p.t_rot-residual_filtration_duration-washing_output.washing_duration)
        end
        fprintf('\nPosition 4\n')
        fprintf('Deliquoring: %2.1f seconds\n',p.t_rot)
    elseif residual_time_pos3 == 0 % the cake does not undergo deliquoring after washing in position 3
        if residual_filtration_duration == 0
            fprintf('\nPosition 3\n')
            fprintf('Washing: %2.1f seconds\n',washing_output.washing_duration)
        else
            fprintf('Filtration: %2.1f seconds\n',residual_filtration_duration)
            fprintf('Washing: %2.1f seconds\n',washing_output.washing_duration)
        end  
        fprintf('\nPosition 4\n')
        fprintf('Deliquoring: %2.1f seconds\n',p.t_rot)
    elseif p.t_rot-abs(residual_time_pos3)>=0 %the washing solvent hold-up finishes to be filtered in position 4
        if residual_filtration_duration == 0
            fprintf('\nPosition 3\n')
            fprintf('Washing: %2.1f seconds\n',p.t_rot)
        else
            fprintf('\nPosition 3\n')
            fprintf('Filtration: %2.1f seconds\n',washing_start)
            fprintf('Washing: %2.1f seconds\n',washing_output.washing_duration-abs(residual_time_pos3))
        end
        fprintf('\nPosition 4\n')
        fprintf('Washing: %2.1f seconds\n',abs(residual_time_pos3))
        fprintf('Deliquoring: %2.1f seconds\n',p.t_rot-abs(residual_time_pos3))
   end
elseif 2*p.t_rot-residual_filtration_duration>0 % no time for washing in position 3
    fprintf('\nPosition 3\n')
    fprintf('Filtration: %2.1f seconds\n',p.t_rot)
    fprintf('\nPosition 4\n')
    fprintf('Filtration: %2.1f seconds\n',washing_start-p.t_rot)
    fprintf('Deliquoring: %2.1f seconds\n',p.t_rot)
else
    fprintf('Filtration not completed within the first 4 positions')
end
fprintf('\nPosition 5\n')
fprintf('Drying: %2.1f seconds\n',p.t_rot)
%% if you want to allow the washing ot go on in position 4, uncomment the following 
% block and lines 58,64-67

% if t_deliq_pos3 > 0
%     fprintf('Washing: %2.1f seconds\n',washing_output.washing_duration)
%     fprintf('Deliquoring: %2.1f seconds\n',t_deliq_pos3)
%     fprintf('\nPosition 4\n')
%     fprintf('Deliquoring: %2.1f seconds\n',p.t_rot)
% else
%     fprintf('Washing: %2.1f seconds\n',washing_output.washing_duration-abs(t_deliq_pos3))
%     fprintf('\nPosition 4\n')
%     fprintf('Deliquoring: %2.1f seconds\n',p.t_rot+abs(t_deliq_pos3))
% end

%% Graphical output
% Filtration variables
t_filt=filt_output.t_filt;
V_filt=filt_output.V_filt;
vol_fr_liq_phase_filt=filt_output.vol_fr_liq_phase_filt;
t_filt_total=filt_output.filtration_duration;

% Pre-deliquoring variables - design charts vs PDE solution

if residual_time_pos2 < 0
    t_pre_deliq=0;
    vol_fr_liq_phase_deliq_pde=filt_output.vol_fr_liq_phase(end);
    vol_fr_liq_phase_charts=filt_output.vol_fr_liq_phase_filt(end);
else
    vol_fr_liq_phase_charts=deliq_output_pos2_charts.vol_fr_liq_phase;
    vol_fr_liq_phase_deliq_pde=output_pos2.vol_fr_liq_phase;
end
solvent_content_vol_eq=p.S_inf*p.E;

% Plot avg solvent content during filtration and deliquoring
figure
plot([t_filt output_pos2.t_deliq+t_filt(end)],[vol_fr_liq_phase_filt vol_fr_liq_phase_charts],...
    [t_filt output_pos2.t_deliq+t_filt(end)],[vol_fr_liq_phase_filt vol_fr_liq_phase_deliq_pde],...
    [0 t_filt(end)+output_pos2.t_deliq(end)],[solvent_content_vol_eq solvent_content_vol_eq])
xlabel('Time [s]')
ylabel('Cake mean vol. solvent content [-]')
set(gca,'fontsize',16,'linewidth',1.3,'xtick',0:50:250)%'xlim',[t_deliq(1) 60],'xtick',0:10:60)
axis([0 p.t_rot 0 vol_fr_liq_phase_filt(end)*1.2] )
lim=get(gca);
lim=lim.YLim;
hold on,plot([t_filt_total t_filt_total],[0 lim(2)],'r','linewidth',.5)
legend('Solvent content - design charts','Solvent content - PDE','Equilibrium moisture content','Beginning deliquoring step')

% Plot solvent content profile after pre-deliquoring
if residual_time_pos2 > 0
    figure
    plot(output_pos2.nodes/p.L_cake,output_pos2.S_final*p.E)
    xlabel('Cake axial coordinate')
    ylabel([{'Cake vol. solvent content profile'},{'@ end of deliquoring [-]'}])
    set(gca,'fontsize',16,'linewidth',1.3,'xtick',0:0.2:1)
end
% Plot solvent content profile after washing
return
figure
plot(washing_output.nodes_washing/p.L_cake,washing_output.xv_mother_liquor)
xlabel('Cake axial coordinate')
ylabel([{'Cake vol. solvent content profile'},{'@ end of washing[-]'}])
set(gca,'fontsize',16,'linewidth',1.3,'xtick',0:0.2:1)


% Plot CSD
% figure
% plot(x,CSD,'linewidth',1.3)
% xlabel('Particle size [m]')
% ylabel('f [#/m^4]')
% set(gca,'fontsize',16,'linewidth',1.3) %,'xlim',[t_deliq(1) 60],'xtick',0:10:60)


