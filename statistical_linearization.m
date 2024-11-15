function [var_displacement, var_velocity, conv, k_eq_time, c_eq_time] =...
    statistical_linearization(m, c, k, M, C, K, freq,...
    time, ndof, epx, q, is_base)

    Mt = M;
    tol = 1e-6;
    maxiter = 30;
    ntime = numel(time);
    nfreq = numel(freq);

    var_displacement = zeros(ndof,ntime);
    var_velocity = zeros(ndof,ntime);
    k_eq_time = zeros(ndof,ntime); 
    c_eq_time = zeros(ndof,ntime);
    keq = zeros(ndof, 1);
    ceq = zeros(ndof, 1);

    if is_base
        M_ps = m'*m;
    else
        M_ps = diag(m);
    end

    for i=1:ntime %loop in time
        t = time(i);
        [Ceq, Keq] = get_equivalent_ck(ceq, keq, ndof);

        sx2 = zeros(ndof, 1);
        sv2 = zeros(ndof, 1);
    
        dkeq = 1000*ones(ndof,1);
        dceq = 1000*ones(ndof,1);
        dkeq_max = max(dkeq);
        dceq_max = max(dceq);

        itera = 0;
        while (dkeq_max > tol || dceq_max > tol) && (itera < maxiter)
            Ct = C + Ceq;
            Kt = K + Keq;
    
            H_ps = zeros(ndof, nfreq);
            H_ps_freq = zeros(ndof, nfreq);
          
            for j=1:nfreq %loop in frequency
                f = freq(j);
                H = get_H(f,Mt,Ct,Kt, q);
                
                ps = evolutionary_power_spectrum(f, t);

                E_ps = ps*M_ps;
                result = real(H*E_ps*H');

                H_ps(:, j) = diag(result);
                H_ps_freq(:, j) = (f.^2)*diag(result);
            end
    
            ceq_1 = ceq;
            keq_1 = keq;
    
            sx2 = zeros(ndof, 1);
            sv2 = zeros(ndof, 1);
            ceq = zeros(ndof, 1);
            keq = zeros(ndof, 1);
         
            for l=1:ndof

                Ex = 2*trapz(freq, H_ps(l,:));% double sided
                Exd = 2*trapz(freq, H_ps_freq(l,:));%double sided
             
                sx2(l) = Ex;
                sv2(l) = Exd;
                ceq(l) = 3*0*c(l)*Exd;
                keq(l) = 3*epx(l)*k(l)*Ex;
 
            end
            
            [Ceq, Keq] = get_equivalent_ck(ceq, keq, ndof);
     
            dceq = abs(ceq - ceq_1)./ceq_1;
            dkeq = abs(keq - keq_1)./keq_1;

            conv(i).ceq(:,itera+1) = abs(ceq - ceq_1);
            conv(i).keq(:,itera+1) = abs(keq - keq_1);

            dceq_max = max(dceq);
            dkeq_max = max(dkeq);
            
            itera = itera + 1;
        end

        var_displacement(:,i) = real(sx2);
        var_velocity(:,i) = real(sv2);
        k_eq_time(:,i) = real(keq);
        c_eq_time(:,i) = real(ceq);
    end

end

function H = get_H(freq, Mt, Ct, Kt, q)
    H = inv(-(freq.^2) * Mt + (1j * freq)^q * Ct + Kt);
end

function [Ceq, Keq]=get_equivalent_ck(ceq, keq, ndof)
    Ceq = zeros(ndof, ndof);
    Keq = zeros(ndof, ndof);

    for i=1:ndof
        Ceq(i, i) = ceq(i);
        Keq(i, i) = keq(i);
    end

end

