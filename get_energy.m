function [energy] = get_energy(varx_sl, omega_eq_2, aux, ndof, time, q)
    for i = 1:numel(ndof)
        for j = 1:numel(omega_eq_2)
            t = time(j);
            sig2t = varx_sl(i,j);
            w2 = omega_eq_2(i,j);
            parfor k = 1:numel(aux)
                Sw = @(x)( evolutionary_power_spectrum(x, t) );
                sfun = (sig2t - 2*integral(@(x)Sw(x)./( abs(w2 - x.^2 + aux(k)*(1i*x).^q).^2 ),0,Inf)).^2;
                energy(j,k) = sfun;
            end
        end
    end
end

