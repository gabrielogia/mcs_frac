function [tt,response_frac] = fde45_mdof_new(fun,tspan,x0,q,wv)
    %N = 10000;
    %tt = linspace(tspan(1), tspan(end), N)';

    tt = tspan;
    N = numel(tspan);
    n = numel(x0);
    ndof =  round(n / 3);
    dt = tt(2) - tt(1);
    
    response_frac = zeros(n,length(tt));
    frac_der = zeros(ndof,length(tt));
    
    response_frac(:,1) = x0;
    x_dot = x0*0;
    response_frac(:,2) = x0 + x_dot*dt;
    
    for i=1:ndof
        dif_x(i,1) = response_frac(i,2)-response_frac(i,1);
    end

    V = length(tt):-1:1;
    U = length(tt)-1:-1:0;
    Tm = V.^(1-q) - U.^(1-q);

    for i=1:length(tt)-1

        w1 = wv(:,i);
        w3 = wv(:,i+1);
        w2 = 0.5*(w1 + w3);
        
        k1 = fun(tt(i),response_frac(:,i),w1);

        soma1 = response_frac(:,i)+dt*k1/2;
        soma1(2*ndof+1:end) = frac_der(:,i);
        k2 = fun(tt(i)+dt/2,soma1,w2);

        soma2 = response_frac(:,i)+dt*k2/2;
        soma2(2*ndof+1:end) = frac_der(:,i);
        k3 = fun(tt(i)+dt/2,soma2,w2);

        soma3 = response_frac(:,i)+dt*k3;
        soma3(2*ndof+1:end) = frac_der(:,i);
        k4 = fun(tt(i)+dt,soma3,w3);
    
        k4(2*ndof+1:end) = frac_der(:,i);
        response_frac(:,i+1) = response_frac(:,i) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
      
        Temp1 = Tm(end-i+1:end);

    
        for ii=1:ndof
            dif_x(ii,i) = response_frac(ii,i+1)-response_frac(ii,i);
            Temp2 = (dif_x(ii,:)*Temp1');
            frac_der(ii,i+1) = (1./(gamma(2-q).*dt.^q)).*Temp2;
            response_frac(2*ndof+ii,i+1) = frac_der(ii,i+1);
        end

    end