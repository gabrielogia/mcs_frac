function P=survival_probability_3(bar,c,time,num_time,beta_eq,weq2,stiff,N, S0, one_integral)
    disp("Computing the survival probability:")
   
    ndof = size(beta_eq,1);
    igmfun = @(et,z)(quadgk(@(t)t.^(et-1).*exp(-t),z,inf,'abstol',1e-12,'reltol',1e-12));
    
    r_2 = zeros(num_time-1, 1);
    P = zeros(ndof, numel(time));
    nr2 = numel(r_2);

    for i=1:N
        cgamma(i) = gamma(i+1);
    end
    
    for dof=1:ndof
        tf=@(z)(interp1(time,pi./sqrt(weq2(dof,:)),z,'pchip'));
        ti=0;

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
        beta_eq_new = interp1(time, beta_eq_dof, time_domain, 'pchip');
        weq2_new = interp1(time, weq2_dof, time_domain, 'pchip');
    
        r_2(1) = 0;
        
        for i = 2:numel(time_domain)
            freq = sqrt(weq2_new(i));
            S = @(t) evolutionary_power_spectrum(freq, t, S0);
            t1 = time_domain(i-1);
            t2 = time_domain(i);
            I1 = integral(@(t) exp(beta_eq_new(i) * t) .* S(t), 0, t1);
            I2 = integral(@(t) exp(beta_eq_new(i) * t) .* S(t), 0, t2);
            r_2(i) = I1/I2;
        end
    
        Q = zeros(nr2,1);
        H = ones(nr2,1);
        F = zeros(nr2,1);
    
        for i = 2:numel(time_domain)
            if(one_integral)
                c_i = c_new(i);
                c_i_1 = c_new(i-1);
                r2 = r_2(i);
                b_i = B/sqrt(c_i);
                b_i_1 = B/sqrt(c_i_1);
                exp1 = exp(-(b_i^2 + b_i_1^2)/(2*(1-r2)));
                I0 = besseli(0,sqrt(r2)*b_i_1*b_i/(1-r2));
                exp2 = exp(-b_i^2/2)*func_int(b_i_1, b_i, r2);
                exp3 = -exp(-b_i_1^2/2)*(1 - func_int(b_i, b_i_1, r2));
                Q(i,1) = exp1*I0 + exp2 + exp3;
                H(i) = 1 - exp((-B^2) / (2 * c_i_1));
                F(i) = Q(i,1) / H(i);

            else
                ci = c_new(i);
                ci1 = c_new(i-1); 
                ri2 = r_2(i);
                A01 = (B^2)/(2*ci*(1-ri2));
                B01 = (B^2)/(2*ci1*(1-ri2));
                D0 = (1-ri2)*exp(-A01)*(1-exp(-B01));
                soma = 0;
    
                for n=1:N
                    AA = (ri2^n); % There is no: (2*n+2)
                    BB = ((ci1*ci)^(n+1))*((1 - ri2)^(2*n+1));
                    CC = prod( (2*(1:n)).^2);% 1:N?
        
                    gi = igmfun(n+1,A01);
                    gi1 = igmfun(n+1,B01);
                    g0 = igmfun(n+1,0);
                    g = cgamma(n);
    
                    Ln = (4^n)*( ci1^(n+1) )*ci*( (1-ri2)^(n+2) )*...
                        ( gi1-g0 )*(-( (ci*(1-ri2))^n )*g +  ( (2*A01)^(-n) )*( g - gi)*B^(2*n) );
         
                    Dn = (AA/(BB*CC))*Ln;
                    
                    soma = soma + Dn;
                end
        
                Q(i) = D0 + soma;
      
                H(i) = 1 - exp((-B^2)/(2*ci1));
                F(i) = Q(i,1)/H(i,1);
            end
        end
     
        F(isnan(F))=0;
    
        Pb_dof = zeros(1, numel(time_domain));
        aux=1;
        for ii=1:numel(time_domain)
            aux = aux*(1 - F(ii));
            Pb_dof(ii) = aux;
        end
    
        P(dof, :) = interp1(time_domain,Pb_dof,time,'pchip');
    end
end

function value = func_int(s1,s2,r2)
    integrand = @(s) (s ./ (1 - r2)) .* ...
        exp(-(s.^2 + r2 * s2^2) / (2 * (1 - r2))) .* ...
        besseli(0, (s .* s2 * sqrt(r2)) / (1 - r2));
    
    value = integral(integrand, 0, s1, 'RelTol', 1e-6, 'AbsTol', 1e-9);
end