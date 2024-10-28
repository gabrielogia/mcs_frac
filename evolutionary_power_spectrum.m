function Sw = evolutionary_power_spectrum(freq, t)
    
    b0 = 0.15;
    S0 = 1;

    aux = (freq./(5*pi)).^2;
    Sw = S0.*aux.*(exp(-b0.*t).*t.^2).*exp(-aux.*t);
    %Sw = S0;
   