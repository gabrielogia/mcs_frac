%% main.m
clc
clear
close all

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

set(groot,'defaultAxesFontSize',12)

%% Structure data
% For reproducibility:
rng(1111);

% Power spectrum (wn, eps)
power_spectrum = "eps";

% Oscillator ('bw', 'duffing')
oscillator = "duffing";

% Number of DOFs:
ndof = 1;

% Fractional derivative:
q = 0.999; 

% Nonlinearity parameter:
epx = 0*ones(1,ndof);

% Mass, damping, and stiffness vectors: 
mass = 1*ones(1,ndof); 
damping = 40*ones(1,ndof);
stiffness = 400*ones(1,ndof);

% Maximum time:
T = 14;

% Barrier:
lam = 0.25;

% Time increment for the Monte Carlo simulation.
dT = 0.0001; %dT = 0.0001;

% Construct matrices M, C, and K:
[M, C, K] = get_mck(mass, damping, stiffness, ndof);

% Maximum frequency of the power spectrum:
fmax_ps = 50; %with bw: 150; 

% Number of samples in the MCS:
ns = 1000;

% Discretization in time and frequency for the Statistical Linearization:
ntime = 300;
nfreq = 1000;

% Run MCS:
run_mcs = true;

% Formulation (optimization, tmdi)
formulation = "optimization"; 

% Base string to save files
str = sprintf('oscillator_%s_ndof_%d_fractional_%.2f_nonlinearity_%.2f_dt_%.4f_mcssamples_%d_damping_%.2f_stiffness_%.2f_formulation_%s_powerspectrum_%s', ...
    oscillator, ndof, q, max(epx), dT, ns, max(damping), max(stiffness), formulation, power_spectrum);

time = linspace(1e-3, T, ntime);

%% Equivalent damping and stiffness
omega_eq_2 = repmat(stiffness,ntime,1)';
beta_eq = repmat(damping,ntime,1)';

%% load montecarlo simulation
load(strcat('data/displacement_variance_', str, '.mat'))

%% Get c(t) by solving the ODE from stochastic averaging.
disp("Solving the ODE to find c(t):")

ic = 0.00000001;
c_integral = zeros(ndof, numel(time));
c_approximation = zeros(ndof, numel(time));
for i=1:ndof
    beta_eq_dof = beta_eq(i,:);
    omega_eq_2_dof = omega_eq_2(i,:);

    beta_eq_dof(1) = findfirstpoint(beta_eq_dof(2),beta_eq_dof(3));
    omega_eq_2_dof(1) = findfirstpoint(omega_eq_2_dof(2),omega_eq_2_dof(3));

    [t, c_integral(i,:)] = ode89(@(t, c_aux) solve_c_mdof(t, c_aux, beta_eq_dof, omega_eq_2_dof, time,q, formulation), time, ic);
    
    formulation = "tmdi";
    [t, c_approximation(i,:)] = ode89(@(t, c_aux) solve_c_mdof(t, c_aux, beta_eq_dof, omega_eq_2_dof, time,q, formulation), time, ic);
end

%% compute analytical c for eps
for i=1:numel(time)
    Sw = @(x)( evolutionary_power_spectrum(x, time(i)) );
    for j=1:ndof
        Sx = @(x)( Sw(x)./( (omega_eq_2(i) - x.^2).^2 + (beta_eq(i)*x).^2 ) );
        s2(i,j) = 2*integral(Sx,0,Inf);
    end
end

%% compute energy diff
k = 2;
[sfun,wq2,bq] = get_energy(s2(end),time(end), omega_eq_2(end), beta_eq(end), q, 101, 101);

idx_w2 = find(round(wq2,5)==round(omega_eq_2(end),5));
idx_beta = find(round(bq,5)==round(beta_eq(end),5));

fig = figure(1);
xlim([0 max(wq2)])
ylim([0 max(bq)])
scatter3(max(stiffness), max(damping), log(sfun(idx_w2, idx_beta)), 100, 'filled')
hold on
surf(wq2, bq, log(sfun))
xlabel('$\omega_{eq}^2(t)$')
ylabel('$\beta_{eq}(t)$')
zlabel('$\Delta$Energy')
colormap jet
shading interp

%% plot omega_eq, beta_eq, and var displacement
fig = figure('color',[1 1 1]);
for i=1:ndof
    subplot(ndof,1,i); 
    hold on
    plot(time, sqrt(s2(:,i)'),'k-','linewidth',2) % SL
    plot(time, sqrt(c_integral(i,:))','r--','linewidth',2) % ODE
    plot(time, sqrt(c_approximation(i,:))','o-','linewidth',1) % ODE
    if run_mcs
        plot(time_out,sqrt(varx_mcs(i,:)),'b:','linewidth',2) % MCS
    end
    legend('Analytical', 'SA - Integral', 'SA - Approximation', 'MCS')
    xlabel('Time (s)')
    ylabel('$\sigma[x(t)] (m)$')
end