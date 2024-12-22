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

% Is Base motion / non-stationary (excitation):
is_base = false;
nonstat = true;

% Number of DOFs:
ndof = 3;

% Power spectrum mag
S0 = 0.2;

% Fractional derivative:
q = 0.50; 

% Nonlinearity parameter:
epx = 1*ones(1,ndof);

% Mass, damping, and stiffness vectors: 
mass = 1*ones(1,ndof); 
damping = 40*ones(1,ndof);
stiffness = 10*damping;

% Bouc-Wen parameters
a_bw = 0.3*ones(1, ndof);
A_bw = 1;
beta_bw = 0.5;
gamma_bw = 0.5;
n_bw = 1;
y0_bw = 0.01;

% Yielding displacement.
xy=y0_bw;

% Maximum time:
T = 4;

% Barrier:
lam = 0.25;

% Time increment for the Monte Carlo simulation.
dT = 1e-3;

% Construct matrices M, C, and K:
if (oscillator == "duffing")
    [M, C, K] = get_mck(mass, damping, stiffness, ndof);
elseif (oscillator == "bw")
    [M, C, K] = get_mck_bw(mass, damping, stiffness, a_bw, ndof, y0_bw);
end

% Maximum frequency of the power spectrum:
fmax_ps = 50; %with bw: 150; 

% Number of samples in the MCS:
ns = 1000;

% Discretization in time and frequency for the Statistical Linearization:
ntime = 200;
nfreq = 1000;

% Formulation (optimization, tmdi)
formulation = "optimization"; 

% Base string to save files
str = sprintf('oscillator_%s_ndof_%d_fractional_%.2f_nonlinearity_%.2f_dt_%.4f_mcssamples_%d_damping_%.2f_stiffness_%.2f_barrier_%.2f_formulation_%s_powerspectrum_%s_S0_%.2f', ...
    oscillator, ndof, q, max(epx), dT, ns, max(damping), max(stiffness), lam, formulation, power_spectrum, S0);

if (oscillator == "bw")
    str = strcat(str, sprintf('_bwparameters_a_%.2f_A_%.2f_beta_%.2f_gamma_%.2f_xy_%.2f', max(a_bw), A_bw, beta_bw, gamma_bw, xy));
end

%% Equivalent damping and stiffness
load(strcat('data/omegaeq_', str, '.mat'));
load(strcat('data/betaeq_', str, '.mat'));

%% load montecarlo simulation
load(strcat('data/displacement_variance_', str, '.mat'))

%% compute energy diff
for i=1:ndof
    [sfun,wq2,bq] = get_energy(varx_sl(i,end), time(end), omega_eq_2(i,end), beta_eq(i,end), q, 71, 71, S0);

    idx_w2 = find(round(wq2,5)==round(omega_eq_2(i,end),5));
    idx_beta = find(round(bq,5)==round(beta_eq(i,end),5));

    data(i).sfun = sfun;
    data(i).wq2 = wq2;
    data(i).bq = bq;
    data(i).idx_w2 = idx_w2;
    data(i).idx_beta = idx_beta;
    
    fig = figure(1);
    subplot(1,ndof,i)
    xlim([0 max(wq2)])
    ylim([0 max(bq)])
    scatter3(omega_eq_2(i,end), beta_eq(i,end), (sfun(idx_w2, idx_beta)), 100, 'filled')
    hold on
    surf(wq2, bq, (sfun))
    xlabel('$\omega_{eq}^2(t)$')
    ylabel('$\beta_{eq}(t)$')
    zlabel('$\Delta$Energy')
    colormap jet
    shading interp
    %view([0,90])
end

saveas(fig, strcat('plots/3ddeltasigma_', str, '.pdf'))

%%
fig = figure(2);
for i=1:ndof
    sfundof = data(i).sfun;
    subplot(ndof,1,i)
    hold on
    plot(bq,(sfundof(data(i).idx_w2, :)))
    scatter(bq(data(i).idx_beta),(sfun(data(i).idx_w2, data(i).idx_beta)),30,'r','filled');
end

saveas(fig, strcat('plots/2ddeltasigma_', str, '.pdf'))



