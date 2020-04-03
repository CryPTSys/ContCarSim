%%%%%%%%%%%%%%%
% reference values for parameteric sensitivities

p.Tg_inlet = 70+273.15; % BC - [K] paracetamol melts at 440 K
p.Vg_N = 100e-6; % 10e-6 for results 1 - flowrate of drying nitrogen in [Nm3/s] - ref: 10NmL/s
p.d_p = 1e-4; % particles mean diameter [m]
p.eps_l_crit = 0.05; % critical solvent content [m3_l/m3]
p.h_M = 500e-10; % 5e-10 for results 1;sensitivity cake height-  solvent mass transfer coefficient [1/s?] - model it!
        % from literature: as long as the driving force is positive, the
        % mass transfer occurs rapidly
% p.h_T_red=1e-5;
p.h_T_value = 100; % W / (m2 K) Bird: 500
p.eps_s = 0.5; % 1 - porosity [m3_s/m3]