clc
close all
clear

%% stiff
load('data/omegaeq_oscillator_soft_ndof_3_fractional_0.50_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.75_powerspectrum_eps_S0_0.20_softparameter_alpha_0.50.mat')

figure(1)
plot(time, omega_eq_2, ':', 'DisplayName','Soft', 'LineWidth',2);
hold on;

load('data/omegaeq_oscillator_duffing_ndof_3_fractional_0.50_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.25_powerspectrum_eps_S0_0.20_duffingparameter_epx_1.00.mat')
plot(time, omega_eq_2, 'DisplayName','Duff', 'LineWidth',2);

load('data/omegaeq_oscillator_bw_ndof_3_fractional_0.50_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.25_powerspectrum_eps_S0_0.20_bwparameters_a_0.50_A_1.00_beta_0.50_gamma_0.50_xy_0.00.mat')
plot(time, omega_eq_2, '--', 'DisplayName','BW', 'LineWidth',2);
legend()
grid on;

%% beta
load('data/betaeq_oscillator_soft_ndof_3_fractional_0.50_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.75_powerspectrum_eps_S0_0.20_softparameter_alpha_0.50.mat')

figure(2)
plot(time, omega_eq_2, ':', 'DisplayName','Soft', 'LineWidth',2);
hold on;

load('data/betaeq_oscillator_duffing_ndof_3_fractional_0.50_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.25_powerspectrum_eps_S0_0.20_duffingparameter_epx_1.00.mat')
plot(time, omega_eq_2, 'DisplayName','Duff', 'LineWidth',2);

load('data/betaeq_oscillator_bw_ndof_3_fractional_0.50_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.25_powerspectrum_eps_S0_0.20_bwparameters_a_0.50_A_1.00_beta_0.50_gamma_0.50_xy_0.00.mat')
plot(time, omega_eq_2, '--', 'DisplayName','BW', 'LineWidth',2);
legend()

grid on