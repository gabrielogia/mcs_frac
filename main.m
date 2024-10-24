%% structure data
clc
clear
close all

%%
ndof = 1;
rng(42);

barrier = 0.02;
is_base = true;
m = 1*ones(1,ndof); %mass vector
damping = 5*ones(1,ndof); %damping vector
stiffness = 100*ones(1,ndof); %stiffness vector
ex = 0.5*ones(1,ndof);
ev = 0*ones(1,ndof);
T = 2;
b0 = 0.15;
S0 = 1;
alpha = 0.75;

[M, C, K] = get_mck(m, damping, stiffness, ndof);

%% plot ps
% Maximum frequency of the excitation process.
fmax_ps = 50; 

vec_time = linspace(0, T, 40);
vec_freq = linspace(0, fmax_ps, 41);

for i=1:numel(vec_time)
    for j=1:numel(vec_freq)
        f = vec_freq(j);
        t = vec_time(i);
        psec(i,j) = evolutionary_power_spectrum(f, t, S0, b0);
    end
end

% figure
% surf(vec_time,vec_freq,psec')

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
ns = 16;
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
[Pb, time_domain] = survival_probability(ndof, barrier, c, beta_eq, omega_eq_2, time, t);

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
    xlim([0, T])
end

str = sprintf('plots/first_passage_time_pdf_ndof_%d_alpha_%.2f_barrier_%.2f_mcs_%d.pdf', ndof, alpha, barrier, ns);
saveas(fig,str)