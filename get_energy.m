function [sfun,weq2,beq] = get_energy(s2t,t,weq2_center,beta_center,q,nw,nb)
    init = 0.01;

    L = init  + rem(weq2_center - init, (weq2_center - init)/nw);
    H = weq2_center*2 - rem(weq2_center*2 - weq2_center, (weq2_center*2 - weq2_center)/nw);
    weq2 = L:(weq2_center - init)/nw:H;

    L = init  + rem(beta_center - init, (beta_center - init)/nb);
    H = beta_center*2 - rem(beta_center*2 - beta_center, (beta_center*2 - beta_center)/nb);
    beq = L:(beta_center - init)/nb:H;

    for i=1:numel(weq2)
        w2 = weq2(i);
        parfor j=1:numel(beq)
            Sw = @(x)( evolutionary_power_spectrum(x, t) );
            sfun(i,j) = (s2t - 2*integral(@(x)Sw(x)./( abs(w2 - x.^2 + beq(j)*(1i*x).^q).^2 ),0,Inf)).^2;
        end
    end
end

