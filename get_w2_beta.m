function [omega_eq_2, beta_eq, beta_original, omega_2_original] = get_w2_beta(ndof, varv_sl, varx_sl, q, dT, T, time, S0, mov)
    omega_eq_2 = varv_sl./varx_sl;

    for i=1:ndof
        omega_eq_2(i,1) = findfirstpoint(omega_eq_2(i,2),omega_eq_2(i,3));
    
        for j=1:numel(time)
            t=time(j);
            sig2t = varx_sl(i,j);

            omega_eq_2_i_j = omega_eq_2(i,j);

            bt = fminbnd(@(y) (sfun(sig2t, omega_eq_2_i_j, t, S0, q, y)), 0.01,1500);
            beq_original(i,j)=bt;

        end
        beq_original(i,1) = findfirstpoint(beq_original(i,2),beq_original(i,3));
        w2 = omega_eq_2(i,:);

        w = sqrt(w2);
        I1 = (time.^(1-q)./(1-q)).*hypergeom(0.5-q/2,[0.5 1.5-q/2],-0.25*w2.*time.^2);
        I2 = (time.^(2-q).*w./(2-q)).*hypergeom(1-q/2,[1.5 2-q/2],-0.25*w2.*time.^2);
        
        if (mov)
            I1 = movmean(I1,30);
            I2 = movmean(I2,30);
        end
        corC = beq_original(i,:).*I1./gamma(1-q);
        corK = beq_original(i,:).*(w.*I2)./gamma(1-q);

        omega_eq_2(i,:) = w2+corK; 
        beq(i,:) = corC;

        beta_original(i,:) = beq_original(i,:);
        omega_2_original(i,:) = w2;
    end
    
    beta_eq = beq;
end