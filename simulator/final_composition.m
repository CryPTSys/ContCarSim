function y = final_composition(p,x,y,n_cycle)
  
    % storing selected states for cakes being processed
    % and calculation of ethanol content in discharged cake (mass fraction)   
    if n_cycle > 0
        if n_cycle>0 && p.stations_working(1)==1
            y.cake_counter(1)=y.cake_counter(1)+1;
            y.states.station1.(['cake_' num2str(y.cake_counter(1))]).t=y.pos1.(['batch_' num2str(n_cycle)]).t;
            y.states.station1.(['cake_' num2str(y.cake_counter(1))]).S=y.pos1.(['batch_' num2str(n_cycle)]).S;
            y.states.station1.(['cake_' num2str(y.cake_counter(1))]).w_EtOH_cake=y.pos1.(['batch_' num2str(n_cycle)]).w_EtOH_cake;

        end
    
        if p.stations_working(2)==1
            y.cake_counter(2)=y.cake_counter(2)+1;
            y.states.station2.(['cake_' num2str(y.cake_counter(2))]).t=y.pos2.(['batch_' num2str(n_cycle-1)]).t;
            y.states.station2.(['cake_' num2str(y.cake_counter(2))]).S=y.pos2.(['batch_' num2str(n_cycle-1)]).S;
            y.states.station2.(['cake_' num2str(y.cake_counter(2))]).w_EtOH_cake=y.pos2.(['batch_' num2str(n_cycle-1)]).w_EtOH_cake;
        end
    
        if p.stations_working(3)==1
            y.cake_counter(3)=y.cake_counter(3)+1;
            y.states.station3.(['cake_' num2str(y.cake_counter(3))]).t=y.pos3.(['batch_' num2str(n_cycle-2)]).t;
            y.states.station3.(['cake_' num2str(y.cake_counter(3))]).S=y.pos3.(['batch_' num2str(n_cycle-2)]).S;
            y.states.station3.(['cake_' num2str(y.cake_counter(3))]).w_EtOH_cake=y.pos3.(['batch_' num2str(n_cycle-2)]).w_EtOH_cake;
        end

        if p.stations_working(4)==1       
            y.cake_counter(4)=y.cake_counter(4)+1;
            y.states.station4.(['cake_' num2str(y.cake_counter(4))]).t=y.pos4.(['batch_' num2str(n_cycle-3)]).t;
            y.states.station4.(['cake_' num2str(y.cake_counter(4))]).S=y.pos4.(['batch_' num2str(n_cycle-3)]).S;
            y.states.station4.(['cake_' num2str(y.cake_counter(4))]).w_EtOH_cake=y.pos4.(['batch_' num2str(n_cycle-3)]).w_EtOH_cake;
            y.states.station4.(['cake_' num2str(y.cake_counter(4))]).Tg_top=y.pos4.(['batch_' num2str(n_cycle-3)]).Tg_top;
            y.final_content(end+1)=...
                y.pos4.(['batch_' num2str(n_cycle-3)]).w_EtOH_cake(end); 
            y.processing_times.(['cake_' num2str(y.cake_counter(4))])=y.sequence.(['batch_' num2str(n_cycle-3)]);

        end
    end    
end