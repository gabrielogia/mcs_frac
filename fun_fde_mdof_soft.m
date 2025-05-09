function x_dot = fun_fde_mdof_soft(t,x,ndof,MiC,MiK,Mi,mass,damping,stiffness,alpha,is_base,w)

    x_dot = zeros(3*ndof,1);
    
    if is_base 
        scale = mass;
    else
        scale = ones(ndof,1);
    end

    H = [zeros(ndof,ndof) eye(ndof) zeros(ndof,ndof);
        -MiK zeros(ndof,ndof) -MiC;
        zeros(ndof,ndof) zeros(ndof,ndof) eye(ndof)];
  
    f = zeros(ndof, 1);
    g0 = zeros(ndof, 1);
    for i=1:ndof
        f(i,1) = scale(i)*w(i);
        g0(i,1) = 0.05*stiffness(i)*abs(x(i))^alpha(i)*sign(x(i));%epx(i)*stiffness(i)*x(i).^3;
    end

    Mf = -Mi*f;
    F = [zeros(ndof,1);Mf;zeros(ndof,1)];
    
    gm = -Mi*g0;
    G = [zeros(ndof,1);gm;zeros(ndof,1)];

    x_dot = H*x + F + G;
