function Sw = evolutionary_power_spectrum(freq, t, S0)
    b0 = 0.2;

    aux1 = (freq./(5*pi)).^2;
    aux2 = (freq./(10*pi)).^2;
    Sw = S0.*aux1.*(exp(-b0.*t).*t.^2).*exp(-aux2.*t);

    % White noise:
    %Sw = S0;