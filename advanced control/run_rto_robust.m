% simulate carousel operation for a given set of decision variables and inputs

function opt_conditions=run_rto_robust(ports_working,x_estim,u,u_nominal,...
    cryst_output_nominal,measurements,n_cycle)
    %% Setup       
    p=rto_carousel_parameters_class;
    p=rto_carousel_parameters(p);
    p.dP=u.P_compr;
    p.Tinlet_drying=u.Tinlet_drying;    
    if ports_working(2)==0
        p.Rm=3e9;
    else
        p.Rm=x_estim.(['batch_' num2str(n_cycle-1)]).Rm+2e9;
    end
    if n_cycle==1
        cryst_output=cryst_output_nominal;
    else
        cryst_output.conc_slurry=measurements.c_slurry_AI101(find(measurements.c_slurry_AI101,1,'last'));      
        cryst_output.T=cryst_output_nominal.T;
    end
    p.T_room=cryst_output.T;
    
    %% Optimization
    options=optimoptions('fmincon','UseParallel',false,'algorithm','sqp','Display','iter',...
        'StepTolerance',1e-11,'OptimalityTolerance',1e-11,'ScaleProblem',true);
    flag=-1;   

    % Pick initial points based on current fouling estimation
    if p.Rm<6e9
        fl0=140e-9;
    elseif p.Rm<10e9
        fl0=110e-9;
    else
        fl0=100e-9;
    end
    
    % Initialization #1
    x0=[fl0, 2e-6];
    while flag <= 0
        [opt_conditions(1,:),~,flag] = fmincon(@obj,x0,[],[],[],[], ...
            [1e-9 1e-6],[5e-7 10e-6],@rto_run_simulation,options,cryst_output,p);
        x0=x0*.9;
    end

    % Initialization #2
    x0=[fl0 4e-6];
    flag=-1;
    while flag <= 0
        [opt_conditions(2,:),~,flag] = fmincon(@obj,x0,[],[],[],[],...
            [1e-9 1e-6],[5e-7 10e-6],@rto_run_simulation,options, cryst_output, p);
        x0=x0*.9;
    end

    % Initialization #3
    x0=[fl0 7e-6];
    flag=-1;
    while flag <= 0
        [opt_conditions(3,:),~,flag] = fmincon(@obj,x0,[],[],[],[],...
            [1e-9 1e-6],[5e-7 10e-6],@rto_run_simulation,options, cryst_output, p);
        x0=x0*.9;
    end
    
    % Initialization #4
    x0=[fl0*0.85 9e-6];
    flag=-1;
    while flag <= 0
        [opt_conditions(4,:),~,flag] = fmincon(@obj,x0,[],[],[],[],...
            [1e-9 1e-6],[5e-7 10e-6],@rto_run_simulation,options, cryst_output, p);
        x0=x0*.9;
    end
    
    % keep the best optimum
    [~,max_index]=max(opt_conditions(:,1));
    opt_conditions=opt_conditions(max_index,:);
    opt_conditions(1)=opt_conditions(2)/opt_conditions(1); 

end

function throughput=obj(x,~,~)
    throughput=-x(1);    
end
