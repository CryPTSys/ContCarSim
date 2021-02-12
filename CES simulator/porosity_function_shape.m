function E = porosity_function_shape(filter_diameter,CSD,x,AR)
    % F. Destro, v1 June 22, 2020
    % Porosity calculation implemented based on Zou and Yu (1996) and Yu, Zou and Standish (1996),
    % Yu and Zou (1998)
    %
    % Inputs:
    % filter_diameter [m]
    % CSD - number based
    % x - [m]
    % kv - volumetric shape coefficient (a dependence on AR and x can be implemented)
    % ks - surface shape coefficient (a dependence on AR and x can be implemented)
    % AR - aspect ratio: considered constant in this version
    kv=0.524;
    ks=3.142;
    sphericity=1;
    
    % bins and grid
    Dx=x(2:end)-x(1:end-1);
    x=(x(2:end)+x(1:end-1))/2;
    CSD=(CSD(2:end)+CSD(1:end-1))/2;
        
    % anonymous functions
    psi=@(kv,ks) pi^(1/3)/ks*(6*kv)^(2/3);
    g=@(r) (1-r).^2+0.4*r.*(1-r).^3.7;
    f=@(r) (1-r).^3.3+2.8*r.*(1-r).^2.7;
    
    %% bin volume fraction, sphericity, equivalent diameters
    X=kv*Dx.*CSD'.*x.^3; % volumes per bin
    X=X/sum(X); % volume fractions    
    
    Vp=kv*x.^3;
    sphericity=1;%psi(kv,ks); % bin sphericity 
    D_eq_sphere=(6*Vp/pi).^(1/3);
    D_eq_pack=D_eq_sphere./sphericity^2.785*exp(2.946*(1-sphericity)); % equivalent packing diameter
    
    %% bin porosity
    % equivalent cylinder
    kv_cyl=pi/4/AR^2;
    ks_cyl=pi*(1/(2*AR^2)+1/AR);
    sphericity_cyl=psi(kv_cyl,ks_cyl);
    Ic=abs(sphericity-sphericity_cyl);
    
    % equivalent disk
    kv_disk=pi/4/AR;
    ks_disk=pi*(0.5+1/AR);
    sphericity_disk=psi(kv_disk,ks_disk);
    Id=abs(sphericity-sphericity_disk);
    
    % initial porosity
    m0=trapz(x,CSD);
    m1=trapz(x,CSD.*x); 
    D_mean=m1./m0;
    E_0_Jeschar=0.375+0.34*D_mean./filter_diameter;
    eps0_cyl=exp(sphericity.^5.58*exp(5.89*(1-sphericity))*log(E_0_Jeschar)); % original paper: 0.40 instead than E_0_Jeschar
    eps0_disk=exp(sphericity.^0.60*exp(0.23*(1-sphericity).^0.45)*log(E_0_Jeschar)); % original paper: 0.40 instead than E_0_Jeschar
    initial_porosity=Id./(Ic+Id).*eps0_cyl+Ic./(Ic+Id)*eps0_disk; 
    
    if AR==1
%         initial_porosity=E_0_Jeschar;
    end
    
    V=1./(1-initial_porosity)*ones(length(x),1);
    
    % calculation mixture porosity
    Vt=zeros(1,length(x));
    
    % discretized approach
    % reverse all vectors - follow apparoach in Yu and Zou (1998), you
    % could avoid doing this inversion
    x=flipud(x);
    D_eq_pack=flipud(D_eq_pack);
    X=flipud(X);
    V=flipud(V);
    
    for i = 1:length(x)
        sum_before=0;
        sum_after=0;
        for j = 1:(i-1)
            r=D_eq_pack(i)/D_eq_pack(j);
            g=(1-r)^2+0.4*r*(1-r)^3.7;
            sum_before = sum_before + X(j)*(V(j)-(V(j)-1)*g-V(i));
        end
        for j = i+1:length(x)
            r=D_eq_pack(j)/D_eq_pack(i);
            f=(1-r)^3.3+2.8*r*(1-r)^2.7;
            sum_after = sum_after + X(j)*(V(j)-V(j)*f-V(i));
        end
        Vt(i)=V(i)+sum_before+sum_after;     
    end
    Vt=max(Vt);
    E=1-1/Vt;
    
%     % continuous approach - switch to volumetric CSD
%     % reverse all vectors
%     for i = 1:length(x)
%         x_smaller=x(1:i-1);
%         x_bigger=x(i+1:end);
%         V_smaller=V(1:i-1);
%         V_bigger=V(i+1:end);
%         CSD_smaller=CSD(1:i-1)';
%         CSD_bigger=CSD(i+1:end)';
%         sum_bigger=0;
%         sum_smaller=0;
%         D_eq_pack_smaller=D_eq_pack(1:i-1);
%         D_eq_pack_bigger=D_eq_pack(i+1:end);
%         
%         if i < length(x)-2
%            sum_bigger = trapz(x_bigger,(V_bigger-(V_bigger-1).*...
%            g(D_eq_pack(i)./D_eq_pack_bigger)-V(i)).*CSD_bigger);
%        end
%        if i > 2   
%            sum_smaller = trapz(x_smaller,(V_smaller-V_smaller.*...
%                f(D_eq_pack_smaller./D_eq_pack(i))-V(i)).*CSD_smaller);
%        end
%         Vt(i)=V(i)+sum_smaller+sum_bigger;    
%         
%     end

 
    