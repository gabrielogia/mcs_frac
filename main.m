%% main.m
clc
clear
close all

%% Structure data 

% For reproducibility:
rng(1111);

% Power spectrum (wn, eps)
power_spectrum = "eps";

% Oscillator ('bw', 'duffing', 'soft')
oscillator = "soft";

% Is Base motion / non-stationary (excitation):
is_base = false;
nonstat = true;

% Number of DOFs:
ndof = 3;

% Power spectrum mag
S0 = 0.2;

% Fractional derivative:
q = 0.50; 

% Duffing Nonlinearity parameter:
epx = 1.0*ones(1,ndof);

% Softening Nonlinearity parameter:
alpha = 0.5*ones(1,ndof);

% Mass, damping, and stiffness vectors: 
mass = 1*ones(1,ndof); 
damping = 40*ones(1,ndof);
stiffness = 10*damping;

% Bouc-Wen parameters
a_bw = 0.5*ones(1, ndof);
A_bw = 1;
beta_bw = 0.5;
gamma_bw = 0.5;
n_bw = 1;
y0_bw = 0.001;

% Yielding displacement.
xy=y0_bw;

% Maximum time:
T = 11;

% Barrier:
lam = 0.25;

% Time increment for the Monte Carlo simulation.
dT = 1e-3;

% Construct matrices M, C, and K:
if (oscillator == "duffing" || oscillator == "soft")
    [M, C, K] = get_mck(mass, damping, stiffness, ndof);
elseif (oscillator == "bw")
    [M, C, K] = get_mck_bw(mass, damping, stiffness, a_bw, ndof, y0_bw);
end

% Maximum frequency of the power spectrum:
fmax_ps = 50;

% Number of samples in the MCS:
ns = 14000;

% Discretization in time and frequency for the Statistical Linearization:
ntime = 400;
nfreq = 2000;

% Base string to save files
str = sprintf('oscillator_%s_ndof_%d_fractional_%.2f_dt_%.4f_mcssamples_%d_damping_%.2f_stiffness_%.2f_barrier_%.2f_powerspectrum_%s_S0_%.2f', ...
    oscillator, ndof, q, dT, ns, max(damping), max(stiffness), lam, power_spectrum, S0);

if (oscillator == "bw")
    str = strcat(str, sprintf('_bwparameters_a_%.2f_A_%.2f_beta_%.2f_gamma_%.2f_xy_%.2f', max(a_bw), A_bw, beta_bw, gamma_bw, xy));
elseif (oscillator == "duffing")
    str = strcat(str, sprintf('_duffingparameter_epx_%.2f', max(epx)));
else
    str = strcat(str, sprintf('_softparameter_alpha_%.2f', max(alpha)));
end

% load mcs previous results
if exist(strcat('data/mcs/mcs_', str, '.mat'))
    load(strcat('data/mcs/mcs_', str, '.mat'))
    run_mcs = false;
else
    run_mcs = true;
end

%% Statistical Linearization
disp(["Running Statistical Linearization:" oscillator])

time = linspace(1e-3, T, ntime);
omega_n = sqrt(eig(inv(M)*K));
freq = linspace(0,fmax_ps,nfreq);

if (oscillator == "duffing")    
    [varx_sl, varv_sl, conv, k_eq, c_eq] = ...
    statistical_linearization(mass, damping, stiffness, M, C, K, freq, time, ndof, epx, q, is_base, S0);
elseif (oscillator == "bw")
    [varx_sl, varv_sl, conv, k_eq, c_eq] = ...
    statistical_linearization_bw(M, C, K, time, A_bw, gamma_bw, beta_bw, fmax_ps, nfreq, q, xy, S0);
elseif (oscillator == "soft")
    [varx_sl, varv_sl, conv, k_eq, c_eq] = ...
    statistical_linearization_soft(mass, damping, stiffness, M, C, K, freq, time, ndof, alpha, q, S0);
end

%% Equivalent damping and stiffness
disp("Getting omega and beta.");
[omega_eq_2, beta_eq, beta_original, w2] = get_w2_beta(ndof, varv_sl, varx_sl, q, dT, T, time, S0, false);

%% Get c(t) by solving the ODE from stochastic averaging.
disp("Solving the ODE to find c(t):")

ic = 0.00000001;
c = zeros(ndof, numel(time));

for i=1:ndof
    beta_eq_dof = beta_eq(i,:);
    omega_eq_2_dof = omega_eq_2(i,:);

    beta_eq_dof(1) = findfirstpoint(beta_eq_dof(2),beta_eq_dof(3));
    omega_eq_2_dof(1) = findfirstpoint(omega_eq_2_dof(2),omega_eq_2_dof(3));

    [t, c(i,:)] = ode89(@(t, c_aux) solve_c_mdof(t, c_aux, beta_eq_dof, omega_eq_2_dof, time,q, S0), time, ic);
end

c = abs(c);

for i=1:ndof
    smaxt(i) = max(c(i,:));
end

smaxi = max(smaxt);
for i=1:ndof
    barrier(i) = lam*sqrt(smaxi);
end

%% Monte Carlo Simulation
if 1
    disp('Running MCS:')
    
    tic
    if (oscillator == "duffing")
    [varx_mcs, time_out, first_passage_time,response,velocity,amplitude] = ...
        monte_carlo(ns,M,C,K,epx,q,mass,damping,stiffness,fmax_ps,...
        nonstat, is_base,T,dT, barrier, S0);
    elseif (oscillator == "soft")
    [varx_mcs, time_out, first_passage_time,response,velocity,amplitude] = ...
        monte_carlo_soft(ns,M,C,K,alpha,q,mass,damping,stiffness,fmax_ps,...
        nonstat, is_base,T,dT, barrier, S0);
    elseif (oscillator == "bw")
        [amplitude, time_out, first_passage_time] = ...
        monte_carlo_bw_new(ns,M,C,K,q,fmax_ps,nonstat,is_base, T, dT, barrier, ndof, A_bw, ...
                           gamma_bw, beta_bw, xy, S0, time, omega_eq_2);
    end
    toc

    save(strcat('data/mcs/mcs_', str, '.mat'), "varx_mcs", "amplitude", "time_out", "first_passage_time")
else
    disp('Jumping MCS')
end

%% plot variance
figure;
plot(time, varx_sl, '--','LineWidth', 2)
hold on;
plot(time_out, varx_mcs, ':', 'LineWidth', 2)
plot(t, c', 'LineWidth', 2)
xlabel('Time')
ylabel('Variance')
legend('DOF 1 - SL', 'DOF 2 - SL', 'DOF 3 - SL', ...
    'DOF 1 - MCS', 'DOF 2 - MCS', 'DOF 3 - MCS', ...
    'DOF 1 - MA', 'DOF 2 - MA', 'DOF 3 - MA', 'Location', 'northwest');
xlim([0 max(time)])

%% First passage
disp("Getting first passage")
amplitude = abs(amplitude);

am = shiftdim(amplitude,1);
amax = max(am(:));
av = linspace(0,amax,200);
ntdim=numel(time_out);

for i=1:ndof
    for j=1:ntdim
        cint = interp1(t,c(i,:),time_out(j),'pchip');
        [pp,xx]=ksdensity(am(j,:,i),av);
        pr(:,j,i)=pp';
        pa(:,j,i)= ((av./cint).*exp(-(av.^2)./(2*cint) ))';
    end
end

%% plot pdf surface
disp('plotting pdf surface')

bar = ones(size(time_out))*barrier(1);
ha = ones(size(time_out))*1000;

fig = figure('color',[1 1 1]);
subplot(2,1,1)
hold on
surf(time_out,av,pr(:,:,end));
shading interp
plot3(time_out,bar,ha,'r','linewidth',2)
view([0,90])
clim([0,200])
xlim([0 4])
xlabel('Time','interpreter','latex', 'FontSize', 14)
ylabel('Amplitude','interpreter','latex', 'FontSize', 14)
title('Empirical probability density function', 'Interpreter', 'latex', 'FontSize', 16)

subplot(2,1,2)
hold on
surf(time_out,av,pa(:,:,end));
shading interp
plot3(time_out,bar,ha,'r','linewidth',2)
view([0,90])
clim([0,200])
xlim([0 4])
xlabel('Time','interpreter','latex', 'FontSize', 14)
ylabel('Amplitude','interpreter','latex', 'FontSize', 14)
title('Analytical probability density function', 'Interpreter', 'latex', 'FontSize', 16)

saveas(fig, strcat('plots/pdfs_', str, '.pdf'))
save(strcat('data/pdfs_', str, '.mat'), "time_out", "av", "pr", "pa", "bar", "ha")

%% plot omega_eq and beta_eq
disp('Plotting omega, beta, and displacement')

fig = figure('color',[1 1 1]);
for i=1:ndof
    subplot(ndof,1,i); 
    hold on
    plot(time,omega_eq_2(i,:),'linewidth',2)
    xlim([0 4])
    xlabel('Time','interpreter','latex', 'FontSize', 14)
    ylabel('$\omega^2_{eq}(t)$','interpreter','latex', 'FontSize', 14)
    title('Oscillator equivalent natural frequency', 'Interpreter', 'latex', 'FontSize', 16)
end

saveas(fig, strcat('plots/omegaeq_', str, '.pdf'))
save(strcat('data/omegaeq_', str, '.mat'), "time", "omega_eq_2", "w2")

fig = figure('color',[1 1 1]);
for i=1:ndof
    subplot(ndof,1,i); 
    hold on
    plot(time, beta_eq(i,:),'linewidth',2)
    xlim([0 4])
    xlabel('Time','interpreter','latex', 'FontSize', 14)
    ylabel('$\beta_{eq}(t)$','interpreter','latex', 'FontSize', 14)
    title('Oscillator equivalent damping', 'Interpreter', 'latex', 'FontSize', 16)
end

saveas(fig, strcat('plots/betaeq_', str, '.pdf'))
save(strcat('data/betaeq_', str, '.mat'), "time", "beta_eq", "beta_original")

%% Survival Probability
bt = beta_eq;
tic
P = survival_probability(barrier, c, time, numel(time), bt, omega_eq_2, 15, S0, oscillator, q);
toc

fig = figure('color',[1 1 1]);
for i=1:ndof
    fpt = first_passage_time(:,i);
    fpt = fpt(fpt>0);
    [fpp,tfp]=ksdensity(fpt,'width',0.1,'Function','survivor');
    subplot(ndof,1,ndof-i+1); 
    hold on
    plot(time, P(i,:)','b','linewidth',2);
    plot(tfp, fpp,'r--','linewidth',2);

    fp_time(i,:) = tfp;
    survival_prob_ksd(i,:) = fpp;

    legend('Analytical', 'MCS')

    title(i)
    xlabel('Time')
    ylabel('Survival Propability')
    xlim([0 4])
    ylim([0 1])
end

saveas(fig, strcat('plots/firsttimepassage_', str, '.pdf'))
save(strcat('data/firsttimepassage_', str, '.mat'), "time", "P", "fp_time", "survival_prob_ksd", "smaxi", "barrier", "lam", "c")