function x_dot = fun_ode_mdof(t,x,ndof,MiC,MiK,Mi,mass,damping,stiffness,epx,is_base,w,tt, oscillator, gamma_bw, beta_bw, A_bw)

    x_dot = zeros(2*ndof,1);
    
    if is_base 
        scale = mass;
    else
        scale = ones(ndof,1);
    end

    H = [zeros(ndof,ndof) eye(ndof); -MiK -MiC];
    f = zeros(ndof, 1);
    g0 = zeros(ndof, 1);

    for i=1:ndof
        f(i,1) = scale(i)*interp1(tt,w(i,:),t,'linear');
        if (oscillator == 'duffing')
            g0(i,1) = epx(i)*stiffness(i)*x(i).^3;
        elseif (oscillator == 'bw')
            dx1 = x(n+i);
            z1 = x(2*n+i);
            g0(2*n+i,1) = (A_bw*dx1 - (beta_bw*dx1(1).*abs(z1).^nn - gamma_bw*abs(dx1).*z1.*abs(dx1).^(nn-1)))/1;
        end
    end

    Mf = -Mi*f;
    F = [zeros(ndof,1);Mf];
    
    gm = -Mi*g0;
    G = [zeros(ndof,1);gm];

    x_dot = H*x + F + G;
