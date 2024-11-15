function [vx,time_out,first_passage_time,state,velo,amplitude] = monte_carlo(ns,M,C,K,epx,q,mass,damping,stiffness,fmax_ps,nonstat,is_base, T, dT, bar) 

    Mi = inv(M);
    MiC = Mi*C;
    MiK = Mi*K;
    ndof = size(M,1);
    
    for i=1:ndof
        EPS{i}.fun = @(f,t)( evolutionary_power_spectrum(f, t) );
    end
    
    tt = 0:dT:T;
    nt = numel(tt);
    
    Nrk = 3000;
    state = zeros(ndof,Nrk,ns);
    velo = zeros(ndof,Nrk,ns);
    wexc = zeros(ndof,nt,ns);
    amplitude = zeros(ndof,Nrk,ns);
    time_out = linspace(tt(1), tt(end), Nrk)';
    
    parfor i=1:ns
        wv = simulate_process(tt, ndof, EPS, nonstat, is_base, fmax_ps);
        xini = zeros(3*ndof,1);

        [t1,X] = fde45_mdof_new(@(t,x,w) fun_fde_mdof(t,x,ndof,MiC,MiK,Mi,...
            mass,damping,stiffness,epx,is_base,w), tt, xini, q, wv);

        for j=1:ndof
            xint = interp1(t1,X(j, :),time_out,'pchip'); %just the displacement
            vint = interp1(t1,X(j+ndof, :),time_out,'pchip'); %just the displacement
            state(j, :, i) = xint;
            velo(j, :, i) = vint;
            amplitude(j,:,i) = get_amplitude(xint);
        end
    end

    [~, n_time, ~] = size(state); % z time dimension.

    for i=1:n_time
        for j=1:ndof
            vx(i, j) = var(state(j, i, :));
        end
    end

    vx = vx';

    %% First passage time
    first_passage_time = zeros(ns,ndof);
    for i=1:ns
        for j=1:ndof
            barrier = bar(j);
            sample_path = amplitude(j,:,i);
            
            time_aux = time_out(abs(sample_path) > barrier);

            if numel(time_aux)==0
                first_passage_time(i,j) = 0;
                warning('Check the barrier');
            else
                first_passage_time(i,j) = time_aux(1);
            end
        end
    end
end

function amp = get_amplitude(x)
    [r,c] = size(x);

    if r>c
        x = x';
    end

    npadd = ceil(max(r,c)/5);
    X = padarray(x',npadd,'symmetric')';

    if r>c
        X = X';
    end
  
    A = abs(hilbert(X));
    amp = A(npadd+1:end-npadd);
end