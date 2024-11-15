function Sw = evolutionary_power_spectrum(freq, t)
    
    b0 = 0.15;
    S0 = 1;

    aux1 = (freq./(5*pi)).^2;
    aux2 = (freq./(10*pi)).^2;
    Sw = S0.*aux2.*(exp(-b0.*t).*t.^2).*exp(-aux2.*t);
    %Sw = S0;
   