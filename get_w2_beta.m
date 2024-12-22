function [omega_eq_2, beta_eq, sfun_value] = get_w2_beta(ndof, varv_sl, varx_sl, q, time, S0)
    omega_eq_2 = varv_sl./varx_sl;
    for i=1:ndof
        omega_eq_2(i,1) = findfirstpoint(omega_eq_2(i,2),omega_eq_2(i,3));
    
        [aux, ~] = optimize(time, i, varx_sl, omega_eq_2, q, S0);

        w2 = omega_eq_2(i,:);
        omega_eq_2(i,:) = w2 + aux(1,:).*w2.^(q).*cos(q*pi/2); 
        beta(i,:) = aux(1,:).*omega_eq_2(i,:).^(q-1).*sin(q*pi/2);

        [aux, sfun_value] = optimize(time, i, varx_sl, omega_eq_2(i,:), q, S0);
        beq(i,:) = aux;
    end
    
    beta_eq = beq;
end

function [beq, sfun_value] = optimize(time, i, varx_sl, w2, q, S0)
    for j=1:numel(time)
        t=time(j);
        sig2t = varx_sl(i,j);

        Sw = @(x)( evolutionary_power_spectrum(x, t, S0) );
        Sx = @(x,y)( Sw(x)./( abs(w2(j) - x.^2 + y*(1i*x).^q).^2 )  );
        sfun = @(y) ( log((sig2t - 2*integral(@(x)Sx(x,y),0,Inf)).^2)  );

        [bt,val] = fminunc(sfun,40);

        beq(j)=bt;
        sfun_value(j)=val;
    end

    beq(1) = findfirstpoint(beq(2),beq(3));
end