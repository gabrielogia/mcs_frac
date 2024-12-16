function [omega_eq_2,beta_eq] = get_w2_beta(formulation, ndof, varv_sl, varx_sl, q, dT, T, time)
    if (formulation == "optimization")
        omega_eq_2 = varv_sl./varx_sl;
        for i=1:ndof
            omega_eq_2(i,1) = findfirstpoint(omega_eq_2(i,2),omega_eq_2(i,3));
        
            for j=1:numel(time)
                t=time(j);
                sig2t = varx_sl(i,j);
        
                Sw = @(x)( evolutionary_power_spectrum(x, t) );
                Sx = @(x,y)( Sw(x)./( abs(omega_eq_2(i,j) - x.^2 + y*(1i*x).^q).^2 )  );
                sfun = @(y) ( (sig2t - 2*integral(@(x)Sx(x,y),0,Inf)).^2  );
        
                bt = fminbnd(sfun,0.01,1500);
                beq(i,j)=bt;
        
            end
            beq(i,1) = findfirstpoint(beq(i,2),beq(i,3));
        
            w2 = omega_eq_2(i,:);
            omega_eq_2(i,:) = w2 + beq(i,:).*w2.^(q).*cos(q*pi/2); 
            beq(i,:) = beq(i,:).*w2.^(q-1).*sin(q*pi/2); 
        end
        
        beta_eq = beq;
    else
        tt = 0:dT:T;
        Nrk = 3000;
        time_out = linspace(tt(1), tt(end), Nrk)';

        for i=1:ndof
            cx = interp1(time,varx_sl(i,:),time_out,'pchip');
            cdx = interp1(time,varv_sl(i,:),time_out,'pchip');
            cxdt = gradient(cx,time_out);

            w2 = cdx./cx;
            EPS = evolutionary_power_spectrum(sqrt(w2),time_out);
            beq(i,:) = (1./cx).*(-cxdt + pi*EPS./w2);

            beta_eq(i,:) = interp1(time_out,beq(i,:),time,'pchip');
            omega_eq_2(i,:) = interp1(time_out,w2,time,'pchip');
        end


    end
end

