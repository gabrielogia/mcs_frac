function [sfun,weq2,beq] = get_energy(s2t,t,weq2max,betamax,q,nw,nb)
    weq2 = linspace(0.1,weq2max,nw);
    beq = linspace(0.1,betamax,nb);

    for i=1:nw
        w2 = weq2(i);
        parfor j=1:nb
            Sw = @(x)( evolutionary_power_spectrum(x, t) );
            sfun(i,j) = (s2t - 2*integral(@(x)Sw(x)./( abs(w2 - x.^2 + beq(j)*(1i*x).^q).^2 ),0,Inf)).^2;
        end
    end
end

