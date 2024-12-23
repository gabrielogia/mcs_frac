function c = sfun(sig2t, omega_eq_2_i_j, t, S0, q, y)    
    Sw = @(x)( evolutionary_power_spectrum(x, t, S0) );
    Sx = @(x,y)( Sw(x)./( abs(omega_eq_2_i_j - x.^2 + y*(1i*x).^q).^2 )  );
    c = ( abs(sig2t - 2*integral(@(x)Sx(x,y),0,Inf)).^2  );
end