function [x,y] = model_drying(t,Dt,p,u,x,y,n_cycle,pos)
    
    % Initialization
    xx=x.(['pos' num2str(pos)]);
        
    % if the liquid saturation is not close to S_inf, you have simultaneous
    % mechanical (i.e. deliquoring) and thermal drying. We pretend they
    % occur in series
    initial_liq_mass_fr_vect=xx.liq_mass_fr_vect;
    %rho_liq=p.rho_liquid_phase_from_mass_fr(mean(initial_liq_mass_fr_vect,2));
%     u.dP=u.dP_drying;
    %xx.S_inf=trapz(p.x,0.155*(1+0.031*p.N_cap_CSD(p.x,rho_liq).^(-0.49)).*p.CSD'/p.m0);    % irreducible cake saturation (minimum saturation that can be achieved by displacement of the interstitial liquid by the applied vacuum    
    SR0=(mean(xx.S)-xx.S_inf)/(1-xx.S_inf);
    residual_deliq=0;
    % if the cake is not deliquored  (SR> 20%), calculate the time needed
    % for getting SR = 20% and carry out deliquoring
    if SR0>0.2 %+0.7*p.S_inf 
        %.dP=u.dP_drying;
        p.Rm=0;
        n_switch=min(max(SR0-.2,0)*100,1);   
        thetaP0=((1-SR0)/(1.08*SR0))^(1/.88)*(1-n_switch)+((1-SR0)/(1.46*SR0))^(1/.48)*n_switch;
        residual_deliq_theta=8-thetaP0;
        visc_liq=p.visc_liquid_phase_from_mass_fr(xx.T,mean(initial_liq_mass_fr_vect,2));
        residual_deliq=min(Dt,residual_deliq_theta*((xx.k.*(u.dP_drying))./(xx.E.*visc_liq*(xx.L_cake.^2).*(1-xx.S_inf)))^-1);  % duration of deliquoring in position 5
        Dt=Dt-residual_deliq;
        u.dP=u.dP_drying;
        [x,y]=model_deliquoring_species_grad(t,residual_deliq,p,u,x,y,n_cycle,pos);
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).t_drying=...
            y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).t_filt;
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).wg=...
            zeros(p.number_volatile_components,length(y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).t_filt));
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).Tg=298*...
            ones(1,length(y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).t_filt));
    end
       
    if Dt > 0    % drying occurs
        
%         x.final_liq_mass_fr_vect=mean(drying_output.final_liq_mass_fr_vect,2);
%         drying_output.vol_cont_impurities=drying_output.final_liq_mass_fr_vect*...
%             p.rho_liquid_phase_from_mass_fr(drying_output.final_liq_mass_fr_vect)./...
%             p.rho_liq_components;
        t_drying_in=unique([0 (p.drying_time_step-round(rem(t,p.drying_time_step),6)):p.drying_time_step:Dt Dt]);
%         t_drying=linspace(0,Dt,100);
    
        % Initial conditions
        xx=x.(['pos' num2str(pos)]);
        initial_liq_mass_fr_vect=xx.liq_mass_fr_vect;       
        S0=xx.S;
        if size(initial_liq_mass_fr_vect,2)>1 && abs(initial_liq_mass_fr_vect(1,1)-initial_liq_mass_fr_vect(1,end))/initial_liq_mass_fr_vect(1)>1e-6     
            initial_liq_mass_fr_vect_rescaled=zeros(p.number_components,xx.number_nodes_drying);
            for i = 1 : size(initial_liq_mass_fr_vect,1)
                initial_liq_mass_fr_vect_rescaled(i,:) = interp1(xx.nodes,initial_liq_mass_fr_vect(i,:),xx.nodes_drying,'pchip');
            end
        else
            initial_liq_mass_fr_vect_rescaled = ones(size(initial_liq_mass_fr_vect,1),...
                xx.number_nodes_drying).*initial_liq_mass_fr_vect(:,1);
        end
        initial_liq_conc=p.conc_from_mass_fr(initial_liq_mass_fr_vect_rescaled);
        if size(S0,2)>1 && abs(S0(1,1)-S0(1,end))/S0(1)>1e-5 
            S0 = interp1(xx.nodes,S0,xx.nodes_drying);
        else
            S0 = ones(1,xx.number_nodes_drying).*S0(:,1);
        end

        % Pre-calculations
        p.dPgdz=u.dP_drying/xx.L_cake;
        p.Pprofile=linspace(101325+u.dP_drying-p.dPgdz*xx.step_grid_drying/2,...
            101325-p.dPgdz*xx.step_grid_drying/2,xx.number_nodes_drying);
        %p.cp_l=mean(p.cp_liquid_phase_from_mass_fr(initial_liq_mass_fr_vect_rescaled));
        %p.rho_l=mean(p.rho_liquid_phase_from_mass_fr(initial_liq_mass_fr_vect_rescaled));
        ug0=u.dP_drying/(p.visc_gas_phase*(xx.alpha*p.rho_sol*xx.L_cake*(1-xx.E)+p.Rm)); % temperature and saturation dependencies added in mex file
        
        rho_g=0.5*(101325+101325-u.dP_drying)/8.314/u.Tinlet_drying*p.MW_N2;
        cp_g=sum(p.cp_N2'.*[1, u.Tinlet_drying, u.Tinlet_drying^2, u.Tinlet_drying^3]);
        Re = rho_g*ug0*xx.m1/xx.m0/p.visc_gas_phase/(1-xx.E*mean(S0)); % bird
        p.h_T = 1000*(2.19.*Re.^(-2/3)).*rho_g.*ug0.*cp_g./(cp_g*p.visc_gas_phase/p.k_N2).^(2/3); % Bird, Re < 50, paragraph 14.5
        p.h_M = [15*1e-5/xx.a_V 15*1e-5/xx.a_V 0]'*1; %             %@ Drying - mass transfer coefficient
    
        % Initial conditions
        epsL_0=S0*xx.E;
        Tg_0=xx.Tg;
        Ts_0=xx.Ts;
        
        
        mass_frG_0=xx.gas_mass_fr_vect;    

        vl_0=initial_liq_conc.*repmat(epsL_0,[size(initial_liq_conc,1),1])./p.rho_liq_components;
        vL_0_vol=vl_0(1:p.number_volatile_components,:);
        vL_0_non_vol=vl_0(p.number_volatile_components+1:end,:);
        epsL_non_vol=sum(vL_0_non_vol,1);
        x0=[Tg_0'; mass_frG_0'; reshape(vL_0_vol',size(vL_0_vol,1)*...
            size(vL_0_vol,2),1)]';

        % sparsity matrix - consistent with the FVM upwind differencing scheme
        S1=triu(tril(ones(xx.number_nodes_drying)),-1); % lower diagonal
        S=repmat(S1,[1+2*p.number_volatile_components,1+2*p.number_volatile_components]);
 
        options = odeset('JPattern',S);
         
        [t_drying,output]=ode15s(@drying_mex_1EB,t_drying_in,x0,options,...
            xx.number_nodes_drying,p.number_volatile_components,epsL_non_vol,...
            u.Tinlet_drying,p.cp_N2,p.rho_sol,xx.E,ug0,p.Pprofile,p.MW_N2,...
            xx.step_grid_drying,p.vl_crit,p.vl_eq,p.coeff_antoine',...
            p.MW_components,p.h_M,xx.a_V,p.rho_liq_components,p.latent_heat,...
            p.cp_liq_components,p.cp_s, p.h_T);
              
        if length(t_drying_in)==2
            t_drying=[t_drying(1) t_drying(end)];
            output=[output(1,:); output(end,:)];
        end
        
        t_drying=t_drying+t+residual_deliq;
        Tgas=output(:,1:xx.number_nodes_drying);
        Tsolid=Tgas;%output(:,xx.number_nodes_drying+1:2*xx.number_nodes_drying);
        mass_frG=output(:,xx.number_nodes_drying+1:end-...
                (p.number_volatile_components)*xx.number_nodes_drying);
        vol_cont_impurities=max(reshape(output(end,end-...
                (p.number_volatile_components)*xx.number_nodes_drying+1:end),...
                xx.number_nodes_drying,p.number_volatile_components)',p.vl_eq(1:p.number_volatile_components));
        vol_cont_impurities=[vol_cont_impurities;
            reshape(vL_0_non_vol,xx.number_nodes_drying,p.number_components-...
            p.number_volatile_components)'];
        eps_l=sum(vol_cont_impurities,1);
        final_liq_conc_vector=vol_cont_impurities./eps_l.*...
            p.rho_liq_components;
        final_liq_mass_fr_vect=final_liq_conc_vector./sum(final_liq_conc_vector);
        S_final=eps_l/xx.E;
        
%         mass_frG_1=mass_frG(:,1:xx.number_nodes_drying)./p.MW_components(1).*...
%             p.MW_N2.*p.Pprofile;
%         mass_frG_2=mass_frG(:,xx.number_nodes_drying+1:xx.number_nodes_drying*2)./...
%             p.MW_components(2).*...
%             p.MW_N2.*p.Pprofile;
%         
%         coeff=p.coeff_antoine(2,:);
%         Tg=298+50;
%         Psat_volatiles = (10.^(coeff(:,1)'-coeff(:,2)'./(coeff(:,3)'+(Tg-273.15)))*133.322); % Antoine Equation - saturation pressure [Pa]
%         drying_driving_force=max(Psat_volatiles-mass_frG_2./...
%             p.MW_components(2)'.*p.MW_N2.*p.Pprofile,0);
%   
        %% Collect outputs
        % states      
        x.(['pos' num2str(pos)]).Tg=Tgas(end,:);
        x.(['pos' num2str(pos)]).Ts=Tsolid(end,:);
        x.(['pos' num2str(pos)]).S=S_final(end,:);
        x.(['pos' num2str(pos)]).nodes=xx.nodes_deliq;
        x.(['pos' num2str(pos)]).liq_mass_fr_vect=final_liq_mass_fr_vect;
        x.(['pos' num2str(pos)]).gas_mass_fr_vect=mass_frG(end,:);

        % measurements
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).t_drying=...
            [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).t_drying,...
            t_drying(2:end)'];
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).Tg=...
            [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).Tg, Tgas(2:end,end)'];
        y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).wg=...
            [y.(['pos' num2str(pos)]).(['cycle_' num2str(n_cycle)]).wg,...
            mass_frG(2:end,[1:p.number_volatile_components]*xx.number_nodes_drying)'];
    end
end
    
 