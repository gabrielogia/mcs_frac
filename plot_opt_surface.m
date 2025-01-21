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
epx = 1.0*ones(1,ndof);

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
y0_bw = 0.01;

% Yielding displacement.
xy=y0_bw;

% Maximum time:
T = 11;

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
ns = 14000;

% Discretization in time and frequency for the Statistical Linearization:
ntime = 200;
nfreq = 1000;

% Base string to save files
str = sprintf('oscillator_%s_ndof_%d_fractional_%.2f_dt_%.4f_mcssamples_%d_damping_%.2f_stiffness_%.2f_barrier_%.2f_powerspectrum_%s_S0_%.2f', ...
    oscillator, ndof, q, dT, ns, max(damping), max(stiffness), lam, power_spectrum, S0);

if (oscillator == "bw")
    str = strcat(str, sprintf('_bwparameters_a_%.2f_A_%.2f_beta_%.2f_gamma_%.2f_xy_%.2f', max(a_bw), A_bw, beta_bw, gamma_bw, xy));
else
    str = strcat(str, sprintf('_duffingparameter_epx_%.2f', max(epx)));
end

%% load simulation results
load(strcat('data/displacement_variance_', str, '.mat'))
load(strcat('data/betaeq_', str, '.mat'))
load(strcat('data/omegaeq_', str, '.mat'))

%% compute surface values
for i=1:ndof
    [sfun_values,wq2,bq] = get_integral_values_mesh(varx_sl(i,end), time(end), w2(i,end), ...
        beta_original(i,end), q, 101, 101, S0);

    idx_w2 = find(round(wq2,5)==round(w2(i,end),5));
    idx_beta = find(round(bq,5)==round(beta_original(i,end),5));

    data(i).sfun_values = sfun_values;
    data(i).wq2 = wq2;
    data(i).bq = bq;
    data(i).idx_w2 = idx_w2;
    data(i).idx_beta = idx_beta;
end    

%% plot
letters = ['a', 'b', 'c'];
for i=1:1
    sfun_values = data(i).sfun_values;
    wq2 = data(i).wq2;
    bq = data(i).bq;
    idx_w2 = data(i).idx_w2;
    idx_beta = data(i).idx_beta;

    fig = figure(1);
    subplot(1,1,i)
    xlim([0 max(wq2)])
    ylim([0 max(bq)])
    zlim([log(min(min(sfun_values))) log(max(max(sfun_values)))])

    scatter3(wq2(idx_w2), bq(idx_beta), log(sfun_values(idx_w2, idx_beta)), 100, 'filled')
    hold on
    surf(wq2, bq, log(sfun_values))
    % aux = sprintf('%s) DOF %d', letters(i), i);
    % title(aux)
    view([18,10])
    %plot3(wq2(idx_w2)*ones(1,1000), linspace(0,max(bq),1000), log(sfun_values(idx_w2, idx_beta))*ones(1,1000), 'r--', 'LineWidth', 2)
    %plot3(linspace(0,max(wq2),1000), bq(idx_beta)*ones(1,1000), log(sfun_values(idx_w2, idx_beta))*ones(1,1000), 'k--', 'LineWidth', 2)

    ax = gca;
    ax.ZMinorTick = 'on';
    ax.ZAxis.TickValues = linspace(log(min(min(sfun_values))), log(max(max(sfun_values))), log(200));

    xlabel('$\omega_{eq}^2(t)$')
    ylabel('$\beta_{eq}(t)$')
    zlabel('$log(\Delta\sigma^2)$')
    colormap jet
    shading interp
end

set(fig,'papersize',[10.5 5.5], 'position', [0 0 800 400]);
print(fig,'plots/3ddeltasigma','-dpdf')

%%
fig = figure(2);
for i=1:1
    sfun_values = data(i).sfun_values;
    bq = data(i).bq;
    subplot(1,1,i)
    hold on
    plot(bq,(sfun_values(data(i).idx_w2, :)))
    scatter(bq(data(i).idx_beta),(sfun_values(data(i).idx_w2, data(i).idx_beta)),50,'r','filled');
    set(gca, 'YScale', 'log')
    % aux = sprintf('%s) DOF %d', letters(i), i);
    % title(aux)
    xlabel('$\beta_{eq}(t)$')
    ylabel('$log(\Delta\sigma^2)$')
    xlim([0, max(bq)])
    grid
end

set(fig,'papersize',[9.5 4.3], 'position', [0 0 750 300]);
print(fig,'plots/2ddeltasigma','-dpdf')