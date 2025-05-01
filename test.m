clc
clear
close all

%%
load('data/firsttimepassage_oscillator_bw_ndof_3_fractional_1.00_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.25_powerspectrum_eps_S0_0.20_bwparameters_a_0.50_A_1.00_beta_0.50_gamma_0.50_xy_0.00.mat')
load('data/omegaeq_oscillator_bw_ndof_3_fractional_1.00_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.25_powerspectrum_eps_S0_0.20_bwparameters_a_0.50_A_1.00_beta_0.50_gamma_0.50_xy_0.00.mat')
load('data/betaeq_oscillator_bw_ndof_3_fractional_1.00_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.25_powerspectrum_eps_S0_0.20_bwparameters_a_0.50_A_1.00_beta_0.50_gamma_0.50_xy_0.00.mat')

%%
% Data
N = 400;
t = linspace(0.01, 11, N)'; % Matches your setup
g = beta_eq(1,:);
s = linspace(0, 5, N)'; % Extended s range

% Rayleigh PDF
K = @(t, s) (s ./ t.^2) .* exp(-s.^2 ./ (2 * t.^2));

% Discretize K(t,s)
K_matrix = zeros(N, N);
for i = 1:N
    K_matrix(i, :) = K(t(i), s)';
end
K_matrix(isinf(K_matrix)) = 1e10; % Cap Infs
K_matrix(isnan(K_matrix)) = 0;

% Integration weights
ds = diff([s; s(end) + (s(end)-s(end-1))]);
K_matrix = K_matrix .* ds';

% Check conditioning
A = K_matrix' * K_matrix;
cond_A = cond(A);
disp(['Condition number: ', num2str(cond_A)]);

% Try lsqnonneg (non-negative f(s))
f = lsqnonneg(K_matrix, g);

% Alternatively, test regularization
lambda = 1e-4; % Adjust based on results
f_reg = (A + lambda * eye(N)) \ (K_matrix' * g);

% Plots
figure;
subplot(3,1,1); plot(t, g, 'b-'); title('g(t)'); xlabel('t'); ylabel('g(t)');

subplot(3,1,2); 
plot(s, f, 'r-', 'DisplayName', 'lsqnonneg'); hold on;
plot(s, f_reg, 'g--', 'DisplayName', 'Regularized'); 
title('f(s)'); xlabel('s'); ylabel('f(s)'); legend;

g_reconstructed = K_matrix * f;
g_reconstructed_reg = K_matrix * f_reg;
subplot(3,1,3); 
plot(t, g, 'b-', 'DisplayName', 'Original'); hold on;
plot(t, g_reconstructed, 'r--', 'DisplayName', 'lsqnonneg');
plot(t, g_reconstructed_reg, 'g:', 'DisplayName', 'Regularized');
legend; title(['Error lsqnonneg: ', num2str(norm(g - g_reconstructed)), ...
    ', Reg: ', num2str(norm(g - g_reconstructed_reg))]); xlabel('t'); ylabel('g(t)');

% Visualize K(t,s)
figure; imagesc(s, t, log10(K_matrix)); colorbar;
xlabel('s'); ylabel('t'); title('log10(K(t,s))');