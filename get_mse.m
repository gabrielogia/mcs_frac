clc
clear
close all

ndof_vec = [3 6 9 12];

for i=1:numel(ndof_vec)
    ndof = ndof_vec(i);
    str = sprintf("oscillator_bw_ndof_%d_fractional_0.75_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.25_powerspectrum_eps_S0_0.20_bwparameters_a_0.50_A_1.00_beta_0.50_gamma_0.50_xy_0.00.mat", ndof);
    load(strcat("data/firsttimepassage_", str))
    
    error = immse(P(:,1:size(survival_prob_ksd,2)), survival_prob_ksd);
    disp(sprintf("NDOF %d; MSE: %.3f", ndof, error));
end