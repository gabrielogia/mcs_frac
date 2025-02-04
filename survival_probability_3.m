function P=survival_probability_3(bar,c,time,num_time,beta_eq,weq2,stiff,N, S0)

disp("Computing the survival probability:")

ndof = size(beta_eq,1);
igmfun = @(et,z)(quadgk(@(t)t.^(et-1).*exp(-t),z,inf,'abstol',1e-12,'reltol',1e-12));

for i=1:N
    cgamma(i) = gamma(i+1);
end

Tf = time(end);
ntf = floor(Tf./(2*pi./(2*sqrt(stiff))));

r_2 = zeros(num_time-1, 1);
Pb = zeros(ndof, num_time);
P = zeros(ndof, numel(time));
nr2 = numel(r_2);

for dof=1:ndof
    tf=@(z)(interp1(time,0.5*2*pi./sqrt(weq2(dof,:)),z,'pchip'));
    ti0 = 0;
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

    c_new = interp1(time, c_dof, time_domain, 'pchip');
    beta_eq_new = interp1(time, beta_eq_dof, time_domain, 'pchip');

    r_2(1) = 0;
    
    for i = 2:numel(time_domain)
        freq = sqrt(weq2(dof,i));
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
        % Constants (example values, change as needed)
        c_i = c_new(i); % example values for c_i and c_{i-1}
        c_i_1 = c_new(i-1);
        r = r_2(i);  % Example value for r
        
        % Laguerre polynomial function
        laguerre_poly = @(n, x) laguerreL(n, x);  % MATLAB built-in Laguerre polynomial
        
        % Define the integrand function
        integrand = @(a_i, a_i_1, c_i, c_i_1, r) ...
            (a_i_1 * a_i / (c_i_1 * c_i)) * ...
            exp(-a_i_1^2 / (2 * c_i_1)) * exp(-a_i^2 / (2 * c_i)) * ...
            sum(arrayfun(@(n) laguerre_poly(n, a_i_1^2 / (2 * c_i_1)) .* ...
            laguerre_poly(n, a_i^2 / (2 * c_i)), 0:N)) * r.^(2 * (0:N));
        
        % Define the bounds for the integration
        lower_bound_a_i_1 = 0;  % Lower limit for a_{i-1}
        upper_bound_a_i_1 = B;  % Upper limit for a_{i-1}
        lower_bound_a_i = B;    % Lower limit for a_i
        upper_bound_a_i = Inf;  % Upper limit for a_i
        
        % Perform the double integral using 'integral2'
        Q(i,1) = integral2(@(a_i_1, a_i) integrand(a_i, a_i_1, c_i, c_i_1, r), ...
            lower_bound_a_i_1, upper_bound_a_i_1, lower_bound_a_i, upper_bound_a_i);

        H(i) = 1 - exp((-B^2)/(2*ci1));
        F(i) = Q(i,1)/H(i,1);
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



