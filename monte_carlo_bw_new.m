function [vx,time_out,first_passage_time,state,velo, hyst, amplitude] = monte_carlo_bw_new(ns,M,C,K,q,fmax_ps,nonstat, ...
                                                                           is_base, T, dT, bar, ndof, A, gamma1, beta1,xy) 
    Ms = M(1:ndof,1:ndof);
    Cs = C(1:ndof,1:ndof);
    Ka = K(1:ndof,1:ndof);
    K1_a = K(1:ndof,ndof+1:2*ndof);

    nn = 1;
    Mi = inv(Ms);
    MiC = Mi*Cs;
    MiKa = Mi*Ka;
    MiK1_a = Mi*K1_a;

    % Excitation Power Spectrum:
    for i=1:ndof
        EPS{i}.fun = @(f,t)( evolutionary_power_spectrum(f, t) );
    end
    
    %MCS
    tt = 0:dT:T;
    nt = numel(tt);
    
    Nrk = 3000;
    state = zeros(ndof,Nrk,ns);
    wexc = zeros(ndof,nt,ns);
    amplitude = zeros(ndof,Nrk,ns);
    time_out = linspace(tt(1), tt(end), Nrk)';
    
    parfor i=1:ns
        wv = simulate_process(tt, ndof, EPS, nonstat, is_base, fmax_ps);
        % Remove it:
        %wv = 2*sin(5*tt);
        xini = zeros(4*ndof,1);
        [t1,X] = fde45_mdof_bw(@(t,x,w) fun_fde_mdof_bw_new(t,x,ndof,MiC,MiKa,MiK1_a,Mi,...
            is_base,w,A,gamma1,beta1,nn, xy), tt, xini, q, wv);

        for j=1:ndof
            xint = interp1(t1,X(j, :),time_out,'pchip');
            vint = interp1(t1,X(j+ndof, :),time_out,'pchip');
            zint = interp1(t1,X(j+2*ndof, :),time_out,'pchip');
            state(j, :, i) = xint;
            velo(j, :, i) = vint;
            hyst(j, :, i) = zint;
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

