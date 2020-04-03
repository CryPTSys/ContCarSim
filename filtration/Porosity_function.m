function [Porosity]=Porosity_function(Filter_diameter,CSD,x)
    % F. Destro, v1 March 27, 2020
    % Porosity calculation implemented based on Ouchiyama, N., Tanaka, T., 1986.Porosity estimation from particle size distribution. Industrial & Engineering Chemistry Fundamentals 25, 125-129.

    x_all=flipud(x(:));
    CSD_all=flipud(CSD(:));
    %     E_0_Dixon=0.4+0.05.*(D_mean./Filter_diameter)+0.412.*(D_mean./Filter_diameter).^2;
    E_denom=zeros(1,length(x_all));    
    % Calculate porosity for the given CSD
    for p = 4:length(x_all)
        xx=x_all(1:p);
        CSD=CSD_all(1:p);
        
        m0=trapz(xx,CSD);
        m1=trapz(xx,CSD.*xx);
        
        D_mean=m1./m0;
        E_0_Jeschar=0.375+0.34*D_mean./Filter_diameter;  % average porosity of packing of uniform sized spheres [-]

        %E_0_deKlerk=0.41+0.35.*exp(0.39.*Filter_diameter./D_mean);

        DD=xx-D_mean;
        DD(DD<=0)=0;

        n_num=sum((xx+D_mean).^2.*(1-3/8*(D_mean./(xx+D_mean))).*CSD);
        n_denom=sum((xx.^3-DD.^3).*CSD);
        n=1+4/13*D_mean*(7-8*E_0_Jeschar)*(n_num./n_denom);
        
        E_denom(p)=sum((DD.^3+(1/n)*((xx+D_mean).^3-DD.^3)).*CSD);
              
    end
    
    E_denom=max(E_denom);
    E_num=sum(x_all.^3.*CSD_all);
    
    Porosity=max(1-E_num/E_denom,0);
