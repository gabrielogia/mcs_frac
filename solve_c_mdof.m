function dc = solve_c_mdof(t, c, beta, omega0, time,q)
    beta_eq = interp1(time, beta, t);

    omega_eq_2 = interp1(time, omega0, t);

    Sw = @(x)( evolutionary_power_spectrum(x, t) );
    
    Sx = @(x)( Sw(x)./( (omega_eq_2 - x.^2).^2 + (beta_eq*x).^2 ) );
    %Sx = @(x)( Sw(x)./( abs(omega_eq_2 - x.^2 + beta_eq*(1i*x).^q).^2 )  );
    s2 = 2*integral(Sx,0,Inf);
    dc = beta_eq.* (s2-c);
    