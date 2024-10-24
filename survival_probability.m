function [Pb, time_domain] = survival_probability(ndof, barrier, c, beta_eq, omega_eq_2, time, t)
    r_2 = zeros(numel(time)-1, 1);
    time_domain = linspace(time(1),time(end),300);
    Pb = zeros(ndof, numel(time_domain));
    
    for dof=1:ndof
        c_dof = c(dof,:);
        beta_eq_dof = beta_eq(dof,:);
        omega_eq_2_dof = omega_eq_2(dof,:);
        beta_eq_dof(1) = beta_eq_dof(2);
        omega_eq_2_dof(1) = omega_eq_2_dof(2);
    
        c_new = interp1(t, c_dof, time_domain, 'pchip');
        beta_eq_new = interp1(t, beta_eq_dof, time_domain, 'pchip');
        omega_eq_2_new = interp1(t, omega_eq_2_dof, time_domain, 'pchip');
        r_2(1) = 0;
        
        for i = 2:numel(time_domain)
            tau = time_domain(i) - time_domain(i-1);
            r_2(i) = (c_new(i-1)/c_new(i))*(1 - (beta_eq_new(i-1))*tau);
        end
        
        N = 5;
        k = 1:1:N;
        
        for i = 2:numel(r_2)
            A = (-barrier^2)/(2*c_new(i)*(1-r_2(i)));
            B = (-barrier^2)/(2*c_new(i-1)*(1-r_2(i)));
            D0 = (1-r_2(i))*exp(A)*(1 - exp(B));
        
            soma = 0;
            for n=1:N
                A = (r_2(i)^n)*((2*n+2));
                B = (c_new(i-1)*c_new(i))^(n+1);
                C = ((1 - r_2(i))^(2*n+1))*prod((2*(1:n)).^2);
                upper_incomplete_gamma_i = igamma(n+1,(barrier^2)/(2*c_new(i)*(1-r_2(i))));
                upper_incomplete_gamma_i_minus = igamma(n+1,(barrier^2)/(2*c_new(i-1)*(1-r_2(i))));
                complete_gamma = gamma(n+1);
                Ln = ((4^n)*((1 - r_2(i))^(2*n+2))*(c_new(i-1)^(n+1))*(c_new(i)^(n+1))*upper_incomplete_gamma_i)*(complete_gamma - upper_incomplete_gamma_i_minus);
                Dn = (A/(B*C))*Ln;
                soma = soma + Dn;
            end
            Q(i,1) = D0 + soma;
            H(i,1) = 1 - exp((-barrier^2)/(2*c_new(i-1)));
            F(i,1) = Q(i,1)/H(i,1);
        end
        
   
        Pb_dof = zeros(1, ndof);
        for i = 1:numel(time_domain)
            aux = 1;
            for p = 1:i
                if (isnan(F(p)))
                    F(p) = 0;
                end
                aux = aux*(1 - F(p));
            end
            Pb_dof(1, i) = aux;
        end
        Pb(dof, :) = Pb_dof;
    end
end