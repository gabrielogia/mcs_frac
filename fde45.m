function [tt,response_frac] = fde45(fun,tspan,x0, q)

N = 15000;
tt = linspace(tspan(1), tspan(end), N)';
n = numel(x0);
% nq = numel(q);
dt = tt(2) - tt(1);

response_frac = zeros(length(tt), n);
frac_der = zeros(length(tt),1);

response_frac(1,:) = x0;
x_dot = fun(0,x0);
response_frac(2,:) = x0 + x_dot*dt;

dif_x(1) = response_frac(2,1)-response_frac(1,1);

%figure
%hold on
for i=2:length(tt)-1
    
    %plot(tt(i),response_frac(i,1),'m.')
    %drawnow

    k1 = fun(tt(i),response_frac(i,:));
    k2 = fun(tt(i)+dt/2,response_frac(i,:)+dt*k1'/2);
    k3 = fun(tt(i)+dt/2,response_frac(i,:)+dt*k2'/2);
    k4 = fun(tt(i)+dt,response_frac(i,:)+dt*k3');

    response_frac(i+1,:) = response_frac(i,:) + (dt/6)*(k1' + 2*k2' + 2*k3' + k4');
    
    for j=1:i-1
    
        Temp1(j) = (i-j).^(1-q) - (i-j-1).^(1-q);
    
    end
    

    dif_x(i-1) = response_frac(i,1)-response_frac(i-1,1);
    
    Temp2 = (dif_x*Temp1');
    frac_der(i) = (1./(gamma(2-q).*dt.^q)).*Temp2;
    response_frac(i+1,end) = frac_der(i);

end




