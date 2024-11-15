function P=survival_probability_3(bar,c,time,num_time,beta_eq,weq2,stiff,N)

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
%figure
for dof=1:ndof


    %time_domain = linspace(time(1),time(end),ntf(dof));
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
        tau = time_domain(i) - time_domain(i-1);
        r_2(i) = (c_new(i-1)/c_new(i))*(1 - (beta_eq_new(i-1))*tau);
    end


    Q = zeros(nr2,1);
    H = ones(nr2,1);
    F = zeros(nr2,1);
    for i = 2:numel(time_domain)
        
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

            %upper_incomplete_gamma_i = igamma(n+1,(barrier^2)/(2*c_new(i)*(1-r_2(i))));
            %upper_incomplete_gamma_i_minus = igamma(n+1,(barrier^2)/(2*c_new(i-1)*(1-r_2(i))));

            gi = igmfun(n+1,A01);
            gi1 = igmfun(n+1,B01);
            g0 = igmfun(n+1,0);
            g = cgamma(n);
           


            Ln = (4^n)*( ci1^(n+1) )*ci*( (1-ri2)^(n+2) )*...
                ( gi1-g0 )*(-( (ci*(1-ri2))^n )*g +  ( (2*A01)^(-n) )*( g - gi)*B^(2*n) );
 
            Li = (B^2)/(2*ci*(1-ri2));

            %Dn = ((-1).^n)*(AA/(BB*CC))*Ln;
            Dn = (AA/(BB*CC))*Ln;
            
            soma = soma + Dn;

        end

        Q(i) = D0 + soma;
  
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



