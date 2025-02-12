function [amplitude, time_out, first_passage_time] = monte_carlo_bw_new(ns,M,C,K,q,fmax_ps,nonstat, ...
                                                                           is_base, T, dT, barrier, ndof, A_bw, ...
                                                                           gamma1, beta1, xy, S0,  ...
                                                                           time, omega_eq_2) 
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
        EPS{i}.fun = @(f,t)( evolutionary_power_spectrum(f, t, S0) );
    end
    
    %MCS
    tt = 0:dT:T;
    Nrk = 3000;
    amplitude = zeros(ndof,Nrk,ns);
    time_out = linspace(tt(1), tt(end), Nrk)';
    tf = zeros(ns,ndof);
    
    parfor i=1:ns
        wv = simulate_process(tt, ndof, EPS, nonstat, is_base, fmax_ps);
        xini = zeros(4*ndof,1);
        [t1,X] = fde45_mdof_bw(@(t,x,w) fun_fde_mdof_bw_new(t,x,ndof,MiC,MiKa,MiK1_a,Mi,...
            is_base,w,A_bw,gamma1,beta1,nn, xy), tt, xini, q, wv);

        for j=1:ndof
            xint = interp1(t1,X(j, :),time_out,'pchip');
            vint = interp1(t1,X(j+ndof, :),time_out,'pchip');

            x = xint;
            dx = vint;
            amp = sqrt(x.^2 + dx.^2./interp1(time,omega_eq_2(j,:),time_out,'pchip'));
            amplitude(j,:,i) = amp;
            time_aux = time_out(abs(amp) > barrier(j));

            if numel(time_aux)==0
                tf(i,j) = NaN;
            else
                tf(i,j) = time_aux(1);
            end
        end
    end

    first_passage_time = tf;
end