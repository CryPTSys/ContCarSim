% Parametric sensitivities

%% results_1
% sens.names = [ {'Tg_inlet'}, {'Vg_N'}, {'d_p'}, {'eps_l_crit'}, {'h_M'},...
%     {'h_T_value'} {'eps_s'} ];
% sens.Tg_inlet = [  50; 60; 70; 80]+273.15;
% sens.Vg_N = [0.5e-6 1e-6 5e-6 10e-6 20e-6 50e-6];
% sens.d_p = [0.5e-5 1e-5 5e-5 1e-4 5e-4];
% sens.eps_l_crit = [0.01 0.03 0.05 0.07 0.1];
% sens.h_M = [0.5e-10 1e-10 5e-10 10e-10 50e-10];
% sens.h_T_value = [0.1 1 10 100 500];
% sens.eps_s = [0.4 0.5 0.6];

% reference values for results_1
% p.Tg_inlet = 70+273.15; % BC - [K] paracetamol melts at 440 K
% p.Vg_N = 10e-6; % flowrate of drying nitrogen in [Nm3/s] - ref: 10NmL/s
% p.d_p = 1e-4; % particles mean diameter [m]
% p.eps_l_crit = 0.05; % critical solvent content [m3_l/m3]
% p.h_M = 5e-10; % solvent mass transfer coefficient [1/s?] - model it!
%         % from literature: as long as the driving force is positive, the
%         % mass transfer occurs rapidly
% % p.h_T_red=1e-5;
% p.Vg_N = 10e-6; % flowrate of drying nitrogen in [Nm3/s] - ref: 10NmL/s
% p.h_T_value = 100; % W / (m2 K) Bird: 500

%% Sensitivity of the impact of cake height on the drying time
% sens.names = [ {'Z'}];
% sens.Z = 1e-3:2.5e-3:0.05;

%% results_3
sens.names = [{'eps_l_crit'},{'h_M'},{'h_T_value'} {'eps_s'} ];
% sens.Tg_inlet = [  70]+273.15;
% sens.Vg_N = [1e-6 50e-6 100e-6 500e-6 1e-5];
% sens.d_p = [0.5e-5 1e-5 5e-5 1e-4 5e-4];
sens.eps_l_crit = linspace(0.01,0.1,10);
sens.h_M = linspace(1e-10,1000e-10,20);
sens.h_T_value = linspace(0.05,100,30);
sens.eps_s = linspace(.4,.5,10);

% Reference values
% p.Tg_inlet = 70+273.15; % BC - [K] paracetamol melts at 440 K
% p.Vg_N = 200e-6; % 10e-6 for results 1 - flowrate of drying nitrogen in [Nm3/s] - ref: 10NmL/s
% p.d_p = 1e-4; % particles mean diameter [m]
% p.eps_l_crit = 0.05; % critical solvent content [m3_l/m3]
% p.h_M = 500e-10; % 5e-10 for results 1;sensitivity cake height-  solvent mass transfer coefficient [1/s?] - model it!
%         % from literature: as long as the driving force is positive, the
%         % mass transfer occurs rapidly
% p.h_T_value = 100; % W / (m2 K) Bird: 500
% p.eps_s = 0.5; % 1 - porosity [m3_s/m3]