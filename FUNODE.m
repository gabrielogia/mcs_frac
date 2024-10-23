function x_dot = FUNODE(t,x,w02,Ca,epx,w,tt)

    x_dot = zeros(3,1); 

    % interpolate excitation.
    f = interp1(tt,w,t,'linear');

    % Duffing with fractional derivative
    x_dot(1) = x(2);
    x_dot(2) = f - (Ca*x(3) + w02*x(1) + epx*w02*x(1).^3);
    x_dot(3) = x(3); % fractional derivative. The last variables will be identified as the fractional derivatives.









