function output=carousel_simulator(cryst_output,p)
        
        % Initialization
        sequence=[]; % vector collecting sequence of processing steps and corresponding duration for the current batch        

        
        %% Position 1: filtration (+ pre-deliquoring)
        % Filtration 
        [filt_output,p]=model_filtration(cryst_output,p);
        sequence{end+1}=['filtration' num2str(filt_output.filtration_duration) ' s'];
         
        % define number of nodes based on cake height
        p.number_nodes_deliq=round(p.L_cake/p.min_length_discr)+1;
        p.number_nodes_filt=p.number_nodes_deliq;
        p.number_nodes_drying=p.number_nodes_deliq;
        p.number_nodes_washing=max(p.number_nodes_deliq,p.number_nodes_washing);

        % Pre-deliquoring
        residual_time_pos1=p.t_rot-filt_output.filtration_duration; % Duration of pre-deliquoring
        if residual_time_pos1 >= 0 % pre-deliquoring occurs in position 1
            output_pos1=model_deliquoring_grad(filt_output,residual_time_pos1,p);   
            residual_filtration_duration=0; % starting time of washing in position 2 
            sequence{end+1}=['pre-deliquoring ' num2str(residual_time_pos1) ' s'];
        elseif residual_time_pos1 < 0 % filtration doesn't finish in position 1 and goes on in position 3
            residual_filtration_duration=abs(residual_time_pos1);
            output_pos1=filt_output;
        end
        
        %% Positions 2-4: (filtration +) washing + deliquoring + drying
        if residual_filtration_duration <= p.t_rot % washing starts in position 2
            washing_output=model_washing(output_pos1,p);
            sequence{end+1}=['washing' num2str(washing_output.washing_duration) ' s'];
            residual_time_pos2=p.t_rot-residual_filtration_duration-washing_output.washing_duration; 
            if residual_time_pos2 >= 0 % the cake undergoes a deliquoring step after washing in position 2
                % Uncomment the following line to see the effect of deliquoring in Position 2 after washing
%                 output_pos3=model_deliquoring_species_grad(washing_output,residual_time_pos2,p); 
                output_pos3=model_deliquoring_species_grad(washing_output,residual_time_pos2+p.t_rot,p);
                sequence{end+1}=['post-deliquoring' num2str(residual_time_pos2+p.t_rot) ' s'];
                drying_duration=p.t_rot;
                drying_output=model_drying(output_pos3,drying_duration,p);
                sequence{end+1}=['drying' num2str(drying_duration) ' s'];
                
            elseif p.t_rot-abs(residual_time_pos2)>=0 %the washing solvent hold-up finishes to be filtered in position 3
                output_pos3=model_deliquoring_species_grad(washing_output,p.t_rot-abs(residual_time_pos2),p);
                sequence{end+1}=['post-deliquoring' num2str(p.t_rot-abs(residual_time_pos2)) ' s'];           
                drying_duration=p.t_rot;
                drying_output=model_drying(output_pos3,drying_duration,p);
                sequence{end+1}=['drying' num2str(drying_duration) ' s'];
                
            elseif 2*p.t_rot-abs(residual_time_pos2)>=0 %the washing solvent hold-up finishes to be filtered in position 4 - no post-deliquoring
                % in this case we have simultaneous thermal drying and deliquoring
                % the drying model by default performs a deliquoring step
                % if the saturation is far from the equilibrium                
                drying_duration=2*p.t_rot-abs(residual_time_pos2);                
                drying_output=model_drying(washing_output,drying_duration,p);
                sequence{end+1}=['drying+deliq' num2str(drying_duration) ' s'];
                
            else % washing does not finish - final composition = final washing comp (not correct, but unfeasible so okay)
                drying_output=washing_output;
                drying_output.final_liq_mass_fr_vect=mean(drying_output.final_liq_mass_fr_vect,2);
                drying_output.vol_cont_impurities=drying_output.final_liq_mass_fr_vect*...
                        p.rho_liquid_phase_from_mass_fr(drying_output.final_liq_mass_fr_vect)./...
                        p.rho_liq_components*p.E;
                drying_output.mass_frG=[];
                drying_output.t_drying=0;
                drying_output.Tgas=[];
                drying_output.Pprofile=[];
                sequence{end+1}='washing not finished';
            end
        elseif residual_filtration_duration<2*p.t_rot % filtration finishes in position3, no time for washing in position 2 
            output_pos3=model_deliquoring_grad(filt_output,p.t_rot-(residual_filtration_duration-p.t_rot),p);
            sequence{end+1}=['post-deliquoring' num2str(p.t_rot-(residual_filtration_duration-p.t_rot)) ' s'];  
            drying_duration=p.t_rot;
            drying_output=model_drying(output_pos3,drying_duration,p);
            sequence{end+1}=['drying' num2str(drying_duration) ' s'];
        elseif residual_filtration_duration<3*p.t_rot % filtration finishes in position 4!
                % in this case we have simultaneous thermal drying and deliquoring
                % for simulation purposes we pretend they occur in series instead
                % than simultaneously
            residual_filtration_duration=residual_filtration_duration-2*p.t_rot; % this residual filtration is to be carried out in position 4
            drying_duration=p.t_rot-residual_filtration_duration; 
            drying_output=model_drying(output_pos1,drying_duration,p);       
            sequence{end+1}=['drying+deliq' num2str(drying_duration) ' s'];
                                             
        else % filtration does not finish
            drying_output.vol_cont_impurities=output_pos1.final_liq_mass_fr_vect*...
                        p.rho_liquid_phase_from_mass_fr(output_pos1.final_liq_mass_fr_vect)./...
                        p.rho_liq_components*p.E;
            sequence{end+1}='filtration not finished';
            drying_output.mass_frG=[];
            drying_output.t_drying=0;
            drying_output.Tgas=[];
            drying_output.Pprofile=[];
        end
        
    %% Create output object
    output.mass_cont_impurities=mean(drying_output.mass_cont_impurities,2);
    output.sequence=sequence;
    output.saturation_before_washing=output_pos1.S_final;
    output.mass_frG=drying_output.mass_frG;
    output.t_drying=drying_output.t_drying;
    output.Tg=drying_output.Tgas;
    output.Pprofile=drying_output.Pprofile;

end
