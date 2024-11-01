function [var_x, var_v, conv, ktime, ctime] =statistical_linearization_bw(M, C, K, time, A, gamma1, beta1, fmax_ps, nfreq)

    tol=1e-6;
    maxiter = 30;

    ndof = numel(M(:,1))/2;
    %Mo = M(1:ndof,1:ndof);
    %Ko = K(1:ndof,1:ndof);
    %wn = sqrt(eig(inv(Mo)*Ko));
    
    freq = linspace(0, fmax_ps, nfreq);
    ntime = numel(time);

    Mt = M;
    
    var_x = zeros(ndof,ntime);
    var_v = zeros(ndof,ntime);
    ktime = zeros(ndof,ntime); 
    ctime = zeros(ndof,ntime); 

    keq = 1+zeros(ndof, 1);
    ceq = 1+1e-101*ones(ndof, 1);

    % Loop in time
    for ii=1:ntime
        t = time(ii);
        
        [~, Ceq, Keq]=create_matrix_bw(ceq, keq, ndof);

        sx2 = zeros(ndof, 1);
        sv2 = zeros(ndof, 1);
    
        dkeq = 1000*ones(ndof,1);
        dceq = 1000*ones(ndof,1);
        dkeq_max = max(dkeq);
        dceq_max = max(dceq);

        itera = 0;
        while (dkeq_max > tol || dceq_max > tol) && itera < maxiter

            Ct = C + Ceq;
            Kt = K + Keq;
    
            H_ps = zeros(2*ndof, nfreq);
            H_ps_freq = zeros(2*ndof, nfreq);
          
            for kk=1:nfreq
                f = freq(kk);
            
                H = getH(f,Mt,Ct,Kt);
               
                aux = zeros(ndof,1);
                auz = zeros(ndof,1);

                for jj = 1:ndof
                    ps = evolutionary_power_spectrum(f, t);
         
                    aux = aux + (2*ps*abs(H(1:ndof, jj)).^2);
                    auz = auz + (2*ps*abs(H(1+ndof:2*ndof, jj)).^2);
                end
                
                H_ps(1:ndof, kk) = aux;
                H_ps_freq(1:ndof, kk) = (f.^2)*aux;
                H_ps(ndof+1:2*ndof,kk) = auz;
                H_ps_freq(ndof+1:2*ndof,kk) = (f.^2)*auz;
            end
    
            ceq_1 = ceq;
            keq_1 = keq;
    
            sx2 = zeros(ndof, 1);
            sv2 = zeros(ndof, 1);
           
            for i=1:ndof
                Ex = trapz(freq, H_ps(i,:));
                Exd = trapz(freq, H_ps_freq(i,:));
                Ez = trapz(freq, H_ps(i + ndof, :));
                Ezd = -(keq(i) / ceq(i)) * Ez;
           
                ceq(i) = (sqrt(2/pi)*(gamma1 * Ezd / sqrt(Exd) + beta1 * sqrt(Ez)) - A);
                keq(i) = (sqrt(2/pi)*(gamma1 * sqrt(Exd) + beta1 * Ezd / sqrt(Ez)));

                sx2(i) = Ex;
                sv2(i) = Exd;
            end
       
            [~, Ceq, Keq]=create_matrix_bw(ceq, keq, ndof);
     
            dceq = abs(ceq - ceq_1)./ceq_1;
            dkeq = abs(keq - keq_1)./keq_1;

            conv(ii).ceq(:,itera+1) = abs(ceq - ceq_1);
            conv(ii).keq(:,itera+1) = abs(keq - keq_1);
     
            dceq_max = max(dceq);
            dkeq_max = max(dkeq);
            
            itera = itera + 1;
        end
        
        var_x(:,ii) = sx2;
        var_v(:,ii) = sv2;
        ktime(:,ii) = keq;
        ctime(:,ii) = ceq;
    end

end

function H = getH(freq, Mt, Ct, Kt)
    H = inv(-(freq.^2) * Mt + 1j * freq * Ct + Kt);
end

function [Meq, Ceq, Keq]=create_matrix_bw(ceq, keq, ndof)
    Meq = zeros(2 * ndof, 2 * ndof);
    Ceq = zeros(2 * ndof, 2 * ndof);
    Keq = zeros(2 * ndof, 2 * ndof);

    cont = 1;
    for i=ndof+1:2 * ndof
        Ceq(i, cont) = ceq(cont);
        Keq(i, i) = keq(cont);
        cont = cont + 1;
    end
end