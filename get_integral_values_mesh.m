function [sfun_values,weq2,beq] = get_integral_values_mesh(s2t,t,weq2_center,beta_center,q,nw,nb,S0)
    k = 0.04;
    L = weq2_center*(1-k);
    H = weq2_center*(1+k);
    weq2 = L:(H - L)/nw:H;
    weq2(numel(weq2)/2) = weq2_center;

    k = 0.07;
    L = beta_center*(1-k);
    H = beta_center*(1+k);
    beq = L:(H - L)/nb:H;
    beq(numel(beq)/2) = beta_center;

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

    sfun_values = zeros(numel(weq2), numel(beq));

    for k=1:numel(weq2)
        w2 = weq2(k);
        parfor j=1:numel(beq)
            beta_eq = beq(j);
            sfun_values(k,j) = sfun(s2t, w2, t, S0, q, beta_eq)
        end
    end
end

