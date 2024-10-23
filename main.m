%% structure data
clc
clear
close all

ndof = 1;
rng(42);

barrier = 0.35;
is_base = true;
m = 1*ones(1,ndof); %mass vector
damping = 1*ones(1,ndof); %damping vector
stiffness = 30*ones(1,ndof); %stiffness vector
ex = 0.5*ones(1,ndof);
ev = 0*ones(1,ndof);
T = 30;
b0 = 0.15;
S0 = 1;
alpha = 0.75;

[M, C, K] = get_mck(m, damping, stiffness, ndof);

%% plot ps
vec_time = linspace(0, T, 40);
vec_freq = linspace(0, 50, 41);

for i=1:numel(vec_time)

    for j=1:numel(vec_freq)

        f = vec_freq(j);
        t = vec_time(i);
        psec(i,j) = evolutionary_power_spectrum(f, t, S0, b0);

    end
end

% figure
% surf(vec_time,vec_freq,psec')

% Maximum frequency of the excitation process.
fmax_ps = 50; 

%% statistical_linearization
disp("sl")
ntime = 100;
nfreq = 1000;
time = linspace(0, T, ntime);
omega_natural = sqrt(eig(inv(M)*K));
freq = linspace(0,max(omega_natural)*3,nfreq);

[var_displacement, var_velocity, conv, k_eq, c_eq] = statistical_linearization(m, damping, stiffness, M, C, K, freq, time, ndof, S0, b0, ex, ev, alpha);

%% stochastic averaging
disp("sa")
omega_eq_2 = var_velocity./var_displacement;
c = var_displacement;

for i=1:ndof
    dc(i,:) = gradient(c(i,:), time);
end
beta_eq = pi.*evolutionary_power_spectrum(sqrt(omega_eq_2), time, S0, b0)./(omega_eq_2.*c) - dc./c;

%% c(t)
disp("c(t)")
epsilon = 0.5;
tspan = time;

ic = 0.00000001;
c_stored = zeros(ndof, numel(time));

for i=1:ndof
    beta_eq_dof = beta_eq(i,:);
    omega_eq_2_dof = omega_eq_2(i,:);

    beta_eq_dof(1) = beta_eq_dof(2);
    omega_eq_2_dof(1) = omega_eq_2_dof(2);

    [t, c_stored(i,:)] = ode45(@(t, c_aux) solve_c_mdof(t, c_aux, beta_eq_dof, omega_eq_2_dof, S0, b0, time), tspan, ic);
end

%% mcs
disp('mcs')
ns = 10;
tic
[var_x, time_out, first_passage_time] = displacement_variance_mcs_mdof(ex, alpha, S0, b0, time, ns, M, C, K, stiffness, ndof, barrier, is_base, fmax_ps);
toc

str = sprintf('data/var_displacements_ndof_%d_alpha_%.2f_barrier_%.2f_mcs_%d.mat', ndof, alpha, barrier, ns);
save(str, 'var_x')
str = sprintf('data/time_out_ndof_%d_alpha_%.2f_barrier_%.2f_mcs_%d.mat', ndof, alpha, barrier, ns);
save(str, 'time_out')
str = sprintf('data/first_passage_time_%d_alpha_%.2f_barrier_%.2f_mcs_%d.mat', ndof, alpha, barrier, ns);
save(str, 'first_passage_time')

%% plot displacements
fig = figure('color',[1 1 1]);
for i=1:ndof
    subplot(ndof,1,i); 
    hold on
    plot(time, c(i,:)','k-','linewidth',2) % SL
    plot(time, c_stored(i,:)','r--','linewidth',2) % ODE
    plot(time_out,var_x(i,:),'b:','linewidth',2) % MCS
    legend('SL', 'SA', 'MCS','interpreter','latex')
    xlabel('Time')
    ylabel('Displacement Variance')
end

str = sprintf('plots/var_displacements_ndof_%d_alpha_%.2f_barrier_%.2f_mcs_%d.pdf', ndof, alpha, barrier, ns);
saveas(fig,str)

%% plot
figure(2)
plot(time,c,time,var_displacement)
title("Oscillator non-stationary response variance E[x^2] = c(t)")
xlabel('Time')
ylabel('Displacement Variance')

figure(3)
plot(t,sqrt(omega_eq_2))
title("\omega_{eq}(t)")
xlabel('Time')
ylabel('Natural frequency equivalent')

figure(4)
plot(t,beta_eq)
title("\beta_{eq}(t)")
xlabel('Time')
ylabel('Damping equivalent')

%% survival probability
disp("sp")

r_2 = zeros(numel(t)-1, 1);
time_domain = linspace(time(1),time(end),300);
%time_domain(1) = time(1);
Pb = zeros(ndof, numel(time_domain));

for dof=1:ndof
    c_dof = c(dof,:);
    beta_eq_dof = beta_eq(dof,:);
    omega_eq_2_dof = omega_eq_2(dof,:);
    beta_eq_dof(1) = beta_eq_dof(2);
    omega_eq_2_dof(1) = omega_eq_2_dof(2);

    %for i = 2:numel(time_domain)
    %    time_domain(i) = time_domain(i-1) + q*2*pi/sqrt((omega_eq_2_dof(i-1)));
    %end

    c_new = interp1(t, c_dof, time_domain, 'pchip');
    beta_eq_new = interp1(t, beta_eq_dof, time_domain, 'pchip');
    omega_eq_2_new = interp1(t, omega_eq_2_dof, time_domain, 'pchip');
    r_2(1) = 0;
    
    for i = 2:numel(time_domain)
        tau = time_domain(i) - time_domain(i-1);
        r_2(i) = (c_new(i-1)/c_new(i))*(1 - (beta_eq_new(i-1))*tau);
    end
    
    N = 5;
    k = 1:1:N;
    
    for i = 2:numel(r_2)
        A = (-barrier^2)/(2*c_new(i)*(1-r_2(i)));
        B = (-barrier^2)/(2*c_new(i-1)*(1-r_2(i)));
        D0 = (1-r_2(i))*exp(A)*(1 - exp(B));
    
        soma = 0;
        for n=1:N
            A = (r_2(i)^n)*((2*n+2));
            B = (c_new(i-1)*c_new(i))^(n+1);
            C = ((1 - r_2(i))^(2*n+1))*prod((2*(1:n)).^2);
            upper_incomplete_gamma_i = igamma(n+1,(barrier^2)/(2*c_new(i)*(1-r_2(i))));
            upper_incomplete_gamma_i_minus = igamma(n+1,(barrier^2)/(2*c_new(i-1)*(1-r_2(i))));
            complete_gamma = gamma(n+1);
            Ln = ((4^n)*((1 - r_2(i))^(2*n+2))*(c_new(i-1)^(n+1))*(c_new(i)^(n+1))*upper_incomplete_gamma_i)*(complete_gamma - upper_incomplete_gamma_i_minus);
            Dn = (A/(B*C))*Ln;
            soma = soma + Dn;
        end
        Q(i,1) = D0 + soma;
        H(i,1) = 1 - exp((-barrier^2)/(2*c_new(i-1)));
        F(i,1) = Q(i,1)/H(i,1);
    end
    
    Pb_dof = zeros(1, ndof);
    for i = 1:numel(time_domain)
        aux = 1;
        for p = 1:i
            if (isnan(F(p)))
                F(p) = 0;
            end
            aux = aux*(1 - F(p));
        end
        Pb_dof(1, i) = aux;
    end
    Pb(dof, :) = Pb_dof;
end

%% plot pb
fig = figure(5);
for i=1:ndof
    [fpp,tfp]=ksdensity(first_passage_time(:,i),'Function','survivor');
    subplot(ndof,1,i); 
    hold on
    plot(time_domain, Pb(i,:)');
    plot(tfp, fpp);
    legend('Analytical', 'Approximation','interpreter','latex')
    title('First passage time')
    xlabel('Time')
    ylabel('Survivor Propability')
    xlim([0, 30])
end

str = sprintf('plots/first_passage_time_pdf_ndof_%d_alpha_%.2f_barrier_%.2f_mcs_%d.pdf', ndof, alpha, barrier, ns);
saveas(fig,str)