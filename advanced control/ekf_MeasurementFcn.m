function yk = ekf_MeasurementFcn(xk, number_nodes)
   
    % Measurement model for drying
    % Inputs:
    %    xk - x[k], states at time k
    %
    % Outputs:
    %    yk - y[k], measurements at time k

%     yk = [xk(number_nodes) xk(number_nodes*3)]; % measurements: outlet temperature
    
    yk = [ xk(number_nodes) ]; % measurements: outlet temperature
%     yk = [ xk(number_nodes) xk(number_nodes*3)]; % measurements: outlet temperature
end