function P=survival_probability(bar,c,time,num_time,beta_eq,N)

disp("Computing the survival probability:")

ndof = size(beta_eq,1);


for i=1:N
    cgamma(i) = gamma(i+1);
end

r_2 = zeros(num_time-1, 1);
Pb = zeros(ndof, num_time);
P = zeros(ndof, numel(time));
nr2 = numel(r_2);
%figure
for dof=1:ndof
    time_domain = linspace(time(1),time(end),num_time);
    
    barrier = bar(dof);
    c_dof = c(dof,:);
    beta_eq_dof = beta_eq(dof,:);

    %beta_eq_dof(1) = beta_eq_dof(2);
    %omega_eq_2_dof = omega_eq_2(dof,:);
    %omega_eq_2_dof(1) = omega_eq_2_dof(2);
    %for i = 2:numel(time_domain)
    %    time_domain(i) = time_domain(i-1) + q*2*pi/sqrt((omega_eq_2_dof(i-1)));
    %end

    c_new = interp1(time, c_dof, time_domain, 'pchip');
    beta_eq_new = interp1(time, beta_eq_dof, time_domain, 'pchip');
    % omega_eq_2_new = interp1(t, omega_eq_2_dof, time_domain, 'pchip');
    r_2(1) = 0;
    
    for i = 2:numel(time_domain)
        tau = time_domain(i) - time_domain(i-1);
        r_2(i) = (c_new(i-1)/c_new(i))*(1 - (beta_eq_new(i-1))*tau);
    end


    Q = zeros(nr2,1);
    H = zeros(nr2,1);
    F = zeros(nr2,1);
    for i = 2:numel(time_domain)

        %A0 = (-barrier^2)/(2*c_new(i)*(1-r_2(i)));
        %B0 = (-barrier^2)/(2*c_new(i-1)*(1-r_2(i)));
        %D0 = (1-r_2(i))*exp(A0)*(1 - exp(B0));

        A01 = -(barrier^2).*(c_new(i-1)+c_new(i))/(2*c_new(i-1)*c_new(i)*(1-r_2(i)));
        B01 = (barrier^2)/(2*c_new(i-1)*(1-r_2(i)));
        D0 = (1-r_2(i))*exp(A01)*(exp(B01)-1);

    
        soma = 0;
        for n=1:N

            A = (r_2(i)^n); % There is no: (2*n+2)
            B = (c_new(i-1)*c_new(i))^(n+1);
            C = ((1 - r_2(i))^(2*n+1))*prod(1*(1:n)).^2;% 1:N?
            
            upper_incomplete_gamma_i = igamma(n+1,(barrier^2)/(2*c_new(i)*(1-r_2(i))));
            upper_incomplete_gamma_i_minus = igamma(n+1,(barrier^2)/(2*c_new(i-1)*(1-r_2(i))));
            complete_gamma = cgamma(n);
         
            Ln = ((4^n)*((1 - r_2(i))^(2*n+2))*...
                (c_new(i-1)^(n+1))*(c_new(i)^(n+1))*...
                upper_incomplete_gamma_i)*...
                (complete_gamma - upper_incomplete_gamma_i_minus);

            % Ln = -c_new(i-1)*c_new(i)*((1 - r_2(i)).^2).*...
            %    ((1./(2*c_new(i-1) - 2*c_new(i)*r_2(i))).^(-n)).*...
            %    ((barrier^2./(2*c_new(i-1) - 2*c_new(i)*r_2(i)) ).^(-n)).*...
            %    (abs(barrier).^(-2*n)).*(complete_gamma - upper_incomplete_gamma_i_minus).*...
            %    (-(barrier.^(2*n)).*((barrier.^2./(c_new(i)*(1 - r_2(i))) ).^n).*...
            %    ((c_new(i)*(1-r_2(i))).^n).*complete_gamma + (barrier^(4*n)) .*(complete_gamma - upper_incomplete_gamma_i) );

            Dn = ((-0.25).^n)*(A/(B*C))*Ln;
            %Dn = (A/(B*C))*Ln;

            soma = soma + Dn;

        end
        Q(i) = D0 + soma;
        H(i) = 1 - exp((-barrier^2)/(2*c_new(i-1)));
        F(i) = Q(i,1)/H(i,1);
    end


    F(isnan(F))=0;
    Pb_dof = zeros(1, numel(time_domain));
    aux=1;
    for ii=1:numel(time_domain)
        
        aux = aux*(1 - F(ii));
        Pb_dof(ii) = aux;

    end
    %for i = 1:numel(time_domain)
    %    aux = 1;
    %    for p = 1:i
    %        if (isnan(F(p)))
    %            F(p) = 0;
    %        end

    %        aux = aux*(1 - F(p) );
            
    %    end

    %    Pb_dof(1, i) = aux;
    %end
    %Pb(dof, :) = Pb_dof;
    P(dof, :) = interp1(time_domain,Pb_dof,time,'pchip');

end



