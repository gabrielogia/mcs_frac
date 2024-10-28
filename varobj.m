function y=varobj(beta_eq,omega_eq_2,time,mass,fmax_ps,nfreq,is_base, ct)

ndof = numel(mass);
freq = linspace(0,fmax_ps,nfreq);

if is_base
    M_ps = mass'*mass;
else
    M_ps = diag(mass);
end

Md = eye(ndof);
Kd = diag(omega_eq_2);
Cd = diag(beta_eq);

for j=1:numel(freq)
    H = inv(-(freq(j).^2) * Md + (1j * freq(j)) * Cd + Kd);
    ps = evolutionary_power_spectrum(freq(j), time);
    E_ps = ps*M_ps;
    result = real(H*E_ps*H');
    H_ps(:, j) = diag(result);
end

for i=1:ndof
    Jvec(i) = abs(ct(i) - 2*trapz(freq,H_ps(i,:))).^2;
end

y=sum(Jvec);
