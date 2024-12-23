function [omega_eq_2, beta_eq, beta_original, omega_2_original] = get_w2_beta(formulation, ndof, varv_sl, varx_sl, q, dT, T, time, S0)
    omega_eq_2 = varv_sl./varx_sl;

    for i=1:ndof
        omega_eq_2(i,1) = findfirstpoint(omega_eq_2(i,2),omega_eq_2(i,3));
    
        for j=1:numel(time)
            t=time(j);
            sig2t = varx_sl(i,j);
    
            % Sw = @(x)( evolutionary_power_spectrum(x, t, S0) );
            % Sx = @(x,y)( Sw(x)./( abs(omega_eq_2(i,j) - x.^2 + y*(1i*x).^q).^2 )  );
            % sfun = @(y) ( abs(sig2t - 2*integral(@(x)Sx(x,y),0,Inf)).^2  );

            omega_eq_2_i_j = omega_eq_2(i,j);

            bt = fminbnd(@(y) (sfun(sig2t, omega_eq_2_i_j, t, S0, q, y)), 0.01,1500);
            beq_original(i,j)=bt;
        end
        beq_original(i,1) = findfirstpoint(beq_original(i,2),beq_original(i,3));
    
        w2 = omega_eq_2(i,:);
        omega_eq_2(i,:) = w2 + beq_original(i,:).*w2.^(q).*cos(q*pi/2); 
        beq(i,:) = beq_original(i,:).*w2.^(q-1).*sin(q*pi/2);

        beta_original(i,:) = beq_original(i,:);
        omega_2_original(i,:) = w2;
    end
    
    beta_eq = beq;
end