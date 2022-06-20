% simulate carousel operation for a given set of decision variables and inputs
function [c, ceq] =rto_run_simulation(x,cryst_output,p)
% Input/Output carousel model, calculating final cake composition for given
% set of inputs

    %% Simulation  
    p.t_rot=x(2)/x(1);
    p.V_slurry=x(2);               

    %% Position 1: filtration (+ pre-deliquoring)
    % Filtration 
    [filt_output,p]=rto_model_filtration(cryst_output,p);   

    % define number of nodes based on cake height
    p.number_nodes=10;

    % Pre-deliquoring
    residual_time_pos1=p.t_rot-filt_output.filtration_duration; % Duration of pre-deliquoring
    if residual_time_pos1 >= 0 % pre-deliquoring occurs in position 1
        output_pos1=rto_model_deliquoring(filt_output,residual_time_pos1,p);   
        residual_filtration_duration=0; % starting time of washing in position 2 
        %sequence{end+1}=['pre-deliquoring ' num2str(residual_time_pos1) ' s'];
    elseif residual_time_pos1 < 0 % filtration doesn't finish in position 1 and goes on in position 3
        residual_filtration_duration=abs(residual_time_pos1);
        output_pos1=filt_output;
    end

    %% Positions 2-4: (filtration + washing )+ deliquoring + drying
    if residual_filtration_duration <= p.t_rot % filtration finished in position 1
        residual_time_pos2=p.t_rot-residual_filtration_duration;
        output_pos3=rto_model_deliquoring(output_pos1,p.t_rot+residual_time_pos2,p);
        drying_duration=p.t_rot;
        drying_output=rto_model_drying(output_pos3,drying_duration,p);
    elseif residual_filtration_duration<2*p.t_rot % filtration finished in position 2
        output_pos3=rto_model_deliquoring(filt_output,p.t_rot-(residual_filtration_duration-p.t_rot),p);
        drying_duration=p.t_rot;
        drying_output=rto_model_drying(output_pos3,drying_duration,p);
    elseif residual_filtration_duration<3*p.t_rot % filtration finishes in position 4!
        % in this case we have simultaneous thermal drying and deliquoring
        % for simulation purposes we pretend they occur in series instead
        % than simultaneously
        residual_filtration_duration=residual_filtration_duration-2*p.t_rot; % this residual filtration is to be carried out in position 4
        drying_duration=p.t_rot-residual_filtration_duration; 
        drying_output=rto_model_drying(output_pos1,drying_duration,p);       
    else % filtration does not finish
        drying_output.mass_avg_profile=1;
    end

    %% Create output object
    c=(drying_output.mass_avg_profile(end)-0.005);
    ceq=[];

end