clc
clear
close all

%% q = 0.999; barrier = 0.25 or 0.75

load(['data/pdfs_oscillator_bw_ndof_3_fractional_1.00_nonlinearity_1.00_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_' ...
    'barrier_0.25_formulation_optimization_powerspectrum_eps_S0_0.20_bwparameters_a_0.30_A_1.00_beta_0.50_gamma_0.50_xy_0.01.mat'])

time = [0.5 1.0 1.5 2.0 2.5 3.0 3.5];
ndof = 3;

for dof = 1:ndof
    plot_pdfs_2d(time_out, time, pa, pr, av, dof)
end

function plot_pdfs_2d(time_out, time, pa, pr, av, dof)
    for i=1:numel(time)
        figure(dof);
        subplot(numel(time), 1, i)
        time_idx = find(time_out >= time(i));
        time_idx = time_idx(1);
        
        pr_time = pr(:,time_idx, dof);
        pa_time = pa(:,time_idx, dof);
        plot(av, pr_time, av, pa_time)
        title(time(i))
        sgtitle(dof)
    end
end