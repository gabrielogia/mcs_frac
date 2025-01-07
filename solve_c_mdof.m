function dc = solve_c_mdof(t, c, beta, omega0, time,q,formulation, S0)
    beta_eq = interp1(time, beta, t);

    omega_eq_2 = interp1(time, omega0, t);

    Sw = @(x)( evolutionary_power_spectrum(x, t, S0) );
    
    if (formulation == "optimization")
        Sx = @(x)( Sw(x)./( (omega_eq_2 - x.^2).^2 + (beta_eq*x).^2 ) );
        s2 = 2*integral(Sx,0,Inf);
    else
        s2 = pi*Sw(sqrt(omega_eq_2))/(omega_eq_2*beta_eq);
    end

    dc = beta_eq.* (s2-c);
    