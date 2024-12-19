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

% Oscillator ('bw', 'duffing')
oscillator = "duffing";

% Number of DOFs:
ndof = 3;

% Fractional derivative:
q = 0.50; 

% Nonlinearity parameter:
epx = 1.4*ones(1,ndof);

% Mass, damping, and stiffness vectors: 
mass = 1*ones(1,ndof); 
damping = 20*ones(1,ndof);
stiffness = 200*ones(1,ndof);

% Maximum time:
T = 14;

% Barrier:
lam = 0.25;

% Time increment for the Monte Carlo simulation.
dT = 0.0010; %dT = 0.0001;

% Construct matrices M, C, and K:
[M, C, K] = get_mck(mass, damping, stiffness, ndof);

% Maximum frequency of the power spectrum:
fmax_ps = 50; %with bw: 150; 

% Number of samples in the MCS:
ns = 12000;

% Discretization in time and frequency for the Statistical Linearization:
ntime = 200;
nfreq = 1000;

% Run MCS:
run_mcs = true;

% Base string to save files
str = sprintf('oscillator_%s_ndof_%d_fractional_%.2f_nonlinearity_%.2f_dt_%.4f_mcssamples_%d_damping_%.2f_stiffness_%.2f', ...
    oscillator, ndof, q, max(epx), dT, ns, max(damping), max(stiffness));

time = linspace(1e-3, T, ntime);

%% Equivalent damping and stiffness
load(strcat('data/omegaeq_', str, '.mat'));
load(strcat('data/betaeq_', str, '.mat'));

%% load montecarlo simulation
load(strcat('data/displacement_variance_', str, '.mat'))

%% compute energy diff
k = 2;

for i=1:ndof
    [sfun,wq2,bq] = get_energy(varx_sl(i,end), time(end), omega_eq_2(i,end), beta_eq(i,end), q, 101, 101);

    idx_w2 = find(round(wq2,5)==round(omega_eq_2(i,end),5));
    idx_beta = find(round(bq,5)==round(beta_eq(i,end),5));
    
    fig = figure(3);
    subplot(1,ndof,i)
    xlim([0 max(wq2)])
    ylim([0 max(bq)])
    scatter3(omega_eq_2(i,end), beta_eq(i,end), log(sfun(idx_w2, idx_beta)), 100, 'filled')
    hold on
    surf(wq2, bq, log(sfun))
    xlabel('$\omega_{eq}^2(t)$')
    ylabel('$\beta_{eq}(t)$')
    zlabel('$\Delta$Energy')
    colormap jet
    shading interp
    %view([0,90])
end