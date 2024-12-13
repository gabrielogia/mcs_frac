function x_dot = fun_fde_mdof_bw_new(t,x,n,MiC,MiKa,MiK1_a,Mi,is_base,w,A,gamma1,beta1,nn,xy)

    x_dot = zeros(4*n,1);
    
    f0 = zeros(n,1);
    H0 = zeros(4*n,1);

    for i=1:n
        f0(i,1) = w(i);
        dx1 = x(n+i);
        z1 = x(2*n+i);
        H0(2*n+i,1) = (A*dx1 - (beta1*dx1(1).*abs(z1).^nn + gamma1*abs(dx1).*z1.*abs(z1).^(nn-1)))/xy;
    end
  
    f0M = Mi*f0;
    F = zeros(4*n,1);
    F(n+1:2*n) = f0M;

    G = zeros(4*n, 4*n);
    G(1:n,n+1:2*n) = eye(n);
    G(n+1:2*n,1:n) = -MiKa;
    G(n+1:2*n,3*n+1:4*n) = -MiC;
    G(n+1:2*n,2*n+1:3*n) = -MiK1_a;
    G(3*n+1:4*n,3*n+1:4*n) = eye(n);

    x_dot = G*x + F + H0;
