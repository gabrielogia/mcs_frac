function P=survival_probability(bar, c, time, num_time, beta_eq, weq2, N, S0, oscillator, q)
    disp("Computing the survival probability:")
   
    ndof = size(beta_eq,1);
    
    P = zeros(ndof, numel(time));
    
    for dof=1:ndof
        ti=0;
        tf = @(z)( interp1(time,1.25*pi./sqrt(weq2(dof,:)), z, 'pchip') );

        cont=2;
        time_domain(1) = ti;
        while ti<=time(end)
            ti = ti + tf(ti);
            time_domain(cont) = ti;
            cont=  cont + 1;
        end
    
        B = bar(dof);
        c_dof = c(dof,:);
        beta_eq_dof = beta_eq(dof,:);
        weq2_dof = weq2(dof,:);
    
        c_new = interp1(time, c_dof, time_domain, 'pchip');
        beta_new = @(t) (interp1(time, beta_eq_dof, t, 'pchip'));
        omega_new = @(t) (interp1(time, sqrt(weq2_dof), t, 'pchip'));
        
        r2 = zeros(numel(time_domain), 1);
        aux=1;
        P_dof = ones(1, numel(time_domain));
        r2(1) = 0;
        nt = numel(time_domain);
        Q = zeros(nt,1);
        H = ones(nt,1);
        F = zeros(nt,1);
    
        for i = 2:numel(time_domain) 
            S = @(t) evolutionary_power_spectrum(omega_new(t), t, S0);
            t1 = time_domain(i-1);
            t2 = time_domain(i);
            I1 = integral(@(t) exp(beta_new(t) .* t) .* S(t), 0, t1);
            I2 = integral(@(t) exp(beta_new(t) .* t) .* S(t), 0, t2);

            c_i = c_new(i);
            c_i_1 = c_new(i-1);
            r2(i) = I1/I2;
            r2_aux = r2(i);

            b_i = B/sqrt(c_i);
            b_i_1 = B/sqrt(c_i_1);
            exp1 = exp(-(b_i^2 + b_i_1^2)/(2*(1-r2_aux)));
            I0 = besseli(0,sqrt(r2_aux)*b_i_1*b_i/(1-r2_aux));
            exp2 = exp(-b_i^2/2)*func_int(b_i_1, b_i, r2_aux);
            exp3 = -exp(-b_i_1^2/2)*(1 - func_int(b_i, b_i_1, r2_aux));

            Q(i,1) = exp1*I0 + exp2 + exp3;

            if isnan(Q(i,1))
                Q(i,1) = Q(i-1,1);
            end


            H(i) = 1 - exp((-B^2) / (2 * c_i_1));
            F(i) = Q(i,1)/H(i,1);
            if isnan(F(i))
                F(i)=0;
            end
            aux = aux*(1 - F(i));
            P_dof(i) = aux;
        end
    
        P(dof, :) = interp1(time_domain, P_dof, time, 'pchip');
    end
end

function value = func_int(s1,s2,r2)
    integrand = @(s) (s ./ (1 - r2)) .* ...
        exp(-(s.^2 + r2 * s2^2) / (2 * (1 - r2))) .* ...
        besseli(0, (s .* s2 * sqrt(r2)) / (1 - r2));
    
    value = integral(integrand, 0, s1, 'RelTol', 1e-6, 'AbsTol', 1e-9);
end