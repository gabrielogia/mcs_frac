function [tt,response_frac] = fde45_mdof(fun,tspan,x0,q)
    %N = 10000;
    %tt = linspace(tspan(1), tspan(end), N)';

    tt = tspan;
    N = numel(tspan);
    n = numel(x0);
    ndof =  round(n / 3);
    nq = numel(q);
    dt = tt(2) - tt(1);
    
    response_frac = zeros(n,length(tt));
    frac_der = zeros(ndof,length(tt));
    
    response_frac(:,1) = x0;
    x_dot = fun(0,x0);
    response_frac(:,2) = x0 + x_dot*dt;
    
    for i=1:ndof
        dif_x(i,1) = response_frac(i,2)-response_frac(i,1);
    end

    for i=2:length(tt)-1
        k1 = fun(tt(i),response_frac(:,i));
        k2 = fun(tt(i)+dt/2,response_frac(:,i)+dt*k1/2);
        k3 = fun(tt(i)+dt/2,response_frac(:,i)+dt*k2/2);
        k4 = fun(tt(i)+dt,response_frac(:,i)+dt*k3);
    
        response_frac(:,i+1) = response_frac(:,i) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
      
        for j=1:i-1
            Temp1(j) = (i-j).^(1-q) - (i-j-1).^(1-q);
        end
    
        for ii=1:ndof
            dif_x(ii,i-1) = response_frac(ii,i)-response_frac(ii,i-1);
            Temp2 = (dif_x(ii,:)*Temp1');
            frac_der(ii,i) = (1./(gamma(2-q).*dt.^q)).*Temp2;
            response_frac(2*ndof+ii,i+1) = frac_der(ii,i);
        end
    end