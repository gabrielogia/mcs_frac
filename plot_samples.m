%% main.m
clc
clear
close all

%% Structure data
tic 

% For reproducibility:
rng(1111);

% Oscillator ('bw', 'duffing')
oscillator = "duffing";

% Is Base motion / non-stationary (excitation):
is_base = false;
nonstat = true;

% Number of DOFs:
ndof = 3;

% Fractional derivative:
q = 0.50; 

% Nonlinearity parameter:
epx = 0.2*ones(1,ndof);

% Mass, damping, and stiffness vectors: 
mass = 1*ones(1,ndof); 
damping = 20*ones(1,ndof);
stiffness = 200*ones(1,ndof);

% Bouc-Wen parameters
a_bw = 0.5*ones(1, ndof);
A_bw = 1;
beta_bw = 0.5;
gamma_bw = 0.5;
n_bw = 1;
y0_bw = 1;

% Yielding displacement.
xy=0.1;

% Maximum time:
T = 20;

% Barrier:
lam = 0.7;

% Time increment for the Monte Carlo simulation.
dT = 0.001; %dT = 0.0001;

% Construct matrices M, C, and K:
if (oscillator == "duffing")
    [M, C, K] = get_mck(mass, damping, stiffness, ndof);
elseif (oscillator == "bw")
    [M, C, K] = get_mck_bw(mass, damping, stiffness, a_bw, ndof, y0_bw);
end

% Maximum frequency of the power spectrum:
fmax_ps = 50; 

% Number of samples in the MCS:
ns = 1;

% Discretization in time and frequency for the Statistical Linearization:
ntime = 200;
nfreq = 1000;

% Run MCS:
run_mcs = true;

% Find the Survival Probability (Yes=1, No=0):
run_fps = true;

% Base string to save files
str = sprintf('oscillator_%s_ndof_%d_fractional_%.2f_nonlinearity_%.2f_dt_%.4f_mcssamples_%d_damping_%.2f_stiffness_%.2f', ...
    oscillator, ndof, q, max(epx), dT, ns, max(damping), max(stiffness));

if (oscillator == "bw")
    str = strcat(str, sprintf('_bwparameters_a_%.2f_A_%.2f_beta_%.2f_gamma_%.2f_xy_%.2f', max(a_bw), A_bw, beta_bw, gamma_bw, xy));
end

%% Statistical Linearization
disp(["Running Statistical Linearization:" oscillator])

time = linspace(1e-3, T, ntime);

if (oscillator == "duffing")
    omega_n = sqrt(eig(inv(M)*K));
    freq = linspace(0,fmax_ps,nfreq);
    
    [varx_sl, varv_sl, conv, k_eq, c_eq] = ...
    statistical_linearization(mass, damping, stiffness, M, C, K, freq, time, ndof, epx, q, is_base);
elseif (oscillator == "bw")
    [varx_sl, varv_sl, conv, k_eq, c_eq] = ...
    statistical_linearization_bw(M, C, K, time, A_bw, gamma_bw, beta_bw, fmax_ps, nfreq, q, xy);
end

%% Equivalent damping and stiffness
omega_eq_2 = varv_sl./varx_sl;

for i=1:ndof
    omega_eq_2(i,1) = findfirstpoint(omega_eq_2(i,2),omega_eq_2(i,3));

    for j=1:numel(time)
        t=time(j);
        sig2t = varx_sl(i,j);
        
        Sw = @(x)( evolutionary_power_spectrum(x, t) );
        Sx = @(x,y)( Sw(x)./( abs(omega_eq_2(i,j) - x.^2 + y*(1i*x).^q).^2 )  );
        sfun = @(y) ( (sig2t - 2*integral(@(x)Sx(x,y),0,Inf)).^2  );

        bt = fminbnd(sfun,0.01,1500);
        beq(i,j)=bt;

    end
    beq(i,1) = findfirstpoint(beq(i,2),beq(i,3));

    w2 = omega_eq_2(i,:);
    omega_eq_2(i,:) = w2 + beq(i,:).*w2.^(q).*cos(q*pi/2); 
    beq(i,:) = beq(i,:).*w2.^(q-1).*sin(q*pi/2); 
end

beta_eq = beq;

%% Get c(t) by solving the ODE from stochastic averaging.
disp("Solving the ODE to find c(t):")

ic = 0.00000001;
c = zeros(ndof, numel(time));
for i=1:ndof
    beta_eq_dof = beta_eq(i,:);
    omega_eq_2_dof = omega_eq_2(i,:);

    beta_eq_dof(1) = findfirstpoint(beta_eq_dof(2),beta_eq_dof(3));
    omega_eq_2_dof(1) = findfirstpoint(omega_eq_2_dof(2),omega_eq_2_dof(3));

    [t, c(i,:)] = ode89(@(t, c_aux) solve_c_mdof(t, c_aux, beta_eq_dof, omega_eq_2_dof, time,q), time, ic);
end

for i=1:ndof
    smaxt(i) = max(c(i,:));
end

smaxi = max(smaxt);
for i=1:ndof
    barrier(i) = lam*sqrt(smaxi);
end

%% Monte Carlo Simulation
if run_mcs
    disp('Running MCS:')

    rng(1111);
    
    if (oscillator == "duffing")
    [varx_mcs, time_out, first_passage_time,response,velocity,amplitude] = ...
        monte_carlo(ns,M,C,K,epx,q,mass,damping,stiffness,fmax_ps,...
        nonstat, is_base,T,dT, barrier);
    elseif (oscillator == "bw")
        [varx_mcs, time_out, first_passage_time,response,amplitude] = ...
        monte_carlo_bw_new(ns,M,C,K,q,fmax_ps,nonstat,is_base, T, dT, barrier, ndof, A_bw, gamma_bw, beta_bw, xy);
    end
end

% %%
% [first_passage_time,amplitude] = time_failure(response,velocity,barrier,omega_eq_2,time_out,time);
% 
% am = shiftdim(amplitude,1);
% amax = max(am(:));
% av = linspace(0,amax,20);
% ntdim=numel(time_out);
% 
% for i=1:ndof
%     for j=1:ntdim
%         cint = interp1(t,c(i,:),time_out(j),'pchip');
%         [pp,xx]=ksdensity(am(j,:,i),av);
%         pr(:,j,i)=pp';
%         pa(:,j,i)= ((av./cint).*exp(-(av.^2)./(2*cint) ))';
%     end
% end
% 
% %% plot pdf surface
% bar = ones(size(time_out))*barrier(1);
% ha = ones(size(time_out))*1000;
% 
% fig = figure('color',[1 1 1]);
% hold on
% surf(time_out,av,pa(:,:,1));
% shading interp
% plot3(time_out,bar,ha,'r','linewidth',2)
% view([0,90])
% clim([0,8])
% xlabel('Time','interpreter','latex', 'FontSize', 25)
% ylabel('Amplitude','interpreter','latex', 'FontSize', 25)
% title('Analytical probability density function', 'Interpreter', 'latex', 'FontSize', 16)
% ylim([0 0.35])
% 
% amp = amplitude(1,:,:);
% 
% B = barrier(1);
% idx = 1:numel(amp(1,:,1));
% h = ones(numel(amp(1,:,1)),1)*100000;
% for i=1:size(amplitude,3)
% 
%     A = amp(1,:,i);
%     id0 = idx(A>B);
%     idxs = 1:id0(1);
% 
%     plot3(time_out(idxs),amp(1,idxs,i),h(idxs),'k');
%     scatter3(time_out(idxs(end)),amp(1,idxs(end),i),h(idxs(end)),100,'g','filled');
% 
% 
% end
% 
% saveas(fig, strcat('plots/pdfs_', str, '.pdf'))
% save(strcat('data/pdfs_', str, '.mat'), "time_out", "av", "pr", "pa")

%%
figure('Color',[1 1 1])

subplot(3,1,1)
plot(time_out,response(3,:), 'Color', [0.9290 0.6940 0.1250], 'linewidth',2)
title('DOF 3')
xlabel('Time (s)', 'interpreter','latex', 'FontSize', 14)
ylabel('Displacement [m]','interpreter','latex', 'FontSize', 14)
aux = sprintf("System response: DOF 3");
title(aux, 'FontSize', 14)
ylim([-0.3, 0.3])
grid(1)

subplot(3,1,2)
plot(time_out,response(2,:), 'Color', [0.8500 0.3250 0.0980],'linewidth',2)
title('DOF 2')
xlabel('Time (s)', 'interpreter','latex', 'FontSize', 14)
ylabel('Displacement [m]','interpreter','latex', 'FontSize', 14)
aux = sprintf("System response: DOF 2");
title(aux, 'FontSize', 14)
ylim([-0.3, 0.3])
grid(1)

subplot(3,1,3)
plot(time_out,response(1,:), 'Color', [0 0.4470 0.7410],'linewidth',2)
title('DOF 1')
xlabel('Time (s)', 'interpreter','latex', 'FontSize', 14)
ylabel('Displacement [m]','interpreter','latex', 'FontSize', 14)
aux = sprintf("System response: DOF 1");
title(aux, 'FontSize', 14)
ylim([-0.3, 0.3])
grid(1)