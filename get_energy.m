function [sfun,weq2,beq] = get_energy(s2t,t,weq2_center,beta_center,q,nw,nb,S0)
    init = 0.01;
    k = 2.5;

    L = init  + rem(weq2_center - init, (weq2_center - init)/nw);
    H = weq2_center*k - rem(weq2_center*k - weq2_center, (weq2_center*k - weq2_center)/nw);
    weq2 = L:(weq2_center - init)/nw:H;

    L = init  + rem(beta_center - init, (beta_center - init)/nb);
    H = beta_center*k - rem(beta_center*k - beta_center, (beta_center*k - beta_center)/nb);
    beq = L:(beta_center - init)/nb:H;

    cond = (size(beq) > size(weq2));
    while cond(2)
        beq(end) = [];
        cond = (size(beq) > size(weq2));
    end

    cond = (size(beq) < size(weq2));
    while cond(2)
        weq2(end) = [];
        cond = (size(beq) < size(weq2));
    end

    for k=1:numel(weq2)
        w2 = weq2(k);
        parfor j=1:numel(beq)
            Sw = @(x)( evolutionary_power_spectrum(x, t, S0) );
            sfun(k,j) = log(abs(s2t - 2*integral(@(x)Sw(x)./( abs(w2 - x.^2 + beq(j)*(1i*x).^q).^2 ),0,Inf)).^2);
        end
    end
end

