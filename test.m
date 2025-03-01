clc
clear
close all

load('data/firsttimepassage_oscillator_bw_ndof_3_fractional_1.00_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.25_powerspectrum_eps_S0_0.20_bwparameters_a_0.50_A_1.00_beta_0.50_gamma_0.50_xy_0.00.mat')
load('data/omegaeq_oscillator_bw_ndof_3_fractional_1.00_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.25_powerspectrum_eps_S0_0.20_bwparameters_a_0.50_A_1.00_beta_0.50_gamma_0.50_xy_0.00.mat')
load('data/betaeq_oscillator_bw_ndof_3_fractional_1.00_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.25_powerspectrum_eps_S0_0.20_bwparameters_a_0.50_A_1.00_beta_0.50_gamma_0.50_xy_0.00.mat')

c_dof = c(1,:);
p_1 = zeros(1, numel(c_dof));

a=0.0001;

for i=1:numel(c_dof)
    p_1(i) = (a/c_dof(i))*exp(-a.^2/(c_dof(i)*2));
end


b = beta_eq(1,:);

B = ifft(fft(b)./fft(p_1));

plot(B)