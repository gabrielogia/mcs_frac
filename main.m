%% structure data
clc
clear
close all

%%
% For reproducibility:
rng(1111);

% Is Base motion / non-stationary (excitation):
is_base = true;
nonstat = true;

% Number of DOFs:
ndof = 2;

% Fractional derivative:
q = 0.5; 

% Nonlinearity parameter:
epx = 0.0*ones(1,ndof);

% Mass, damping, and stiffness vectors: 
mass = 1*ones(1,ndof); 
damping = 10*ones(1,ndof);
stiffness = 100*ones(1,ndof);

% Maximum time:
T = 5;

% Barrier:
barrier = 0.01;

% Time increment for the Monte Carlo simulation.
dT = 0.01; %dT = 0.0001;

% Construct matrices M, C, and K:
[M, C, K] = get_mck(mass, damping, stiffness, ndof);

% Maximum frequency of the power spectrum:
fmax_ps = 50; 

% Number of samples in the MCS:
ns = 240;

% Discretization in time and frequency for the Statistical Linearization:
ntime = 100;
nfreq = 1000;

%% Monte Carlo Simulation

run_mcs = input('Run Monte Carlo Simulation (Yes=1, No=0):');

if run_mcs
    disp('Running MCS:')
    
    [varx_mcs, time_out, first_passage_time] = ...
        monte_carlo(ns,M,C,K,epx,q,mass,damping,stiffness,fmax_ps,...
        nonstat, is_base,T,dT, barrier);

end

%% Statistical Linearization
disp("Running Statistical Linearization:")

time = linspace(0, T, ntime);
omega_n = sqrt(eig(inv(M)*K));
freq = linspace(0,fmax_ps,nfreq);

[varx_sl, varv_sl, conv, k_eq, c_eq] = ...
    statistical_linearization(mass, damping, stiffness, M, C, K,...
    freq, time, ndof, epx, q, is_base);

%% Equivalent damping and stiffness

omega_eq_2 = varv_sl./varx_sl;
omega_eq_2(:,1) = omega_eq_2(:,2);

for i=1:ndof

    for j=1:numel(time)
        t=time(j);
        sig2t = varx_sl(i,j);
        Sw = @(x)( evolutionary_power_spectrum(x, t) );
        Sx = @(x,y)( Sw(x)./( (omega_eq_2(i,j) - x.^2).^2 + (y*x).^2 )  );
        sfun = @(y) ( (sig2t - 2*integral(@(x)Sx(x,y),0,Inf)).^2  );

        bt = fminbnd(sfun,0.01,500);
        beq(i,j)=bt;

    end
    beq(i,1) = beq(i,2);

end

beta_eq = beq;


%%
figure('color',[1 1 1]);
for i=1:ndof
    subplot(ndof,1,i); 
    hold on
    plot(time,omega_eq_2(i,:),'linewidth',2) % MCS
    xlabel('Time','interpreter','latex')
    ylabel('$\omega^2_{eq}(t)$','interpreter','latex')
end

figure('color',[1 1 1]);
for i=1:ndof
    subplot(ndof,1,i); 
    hold on
    plot(time, beta_eq(i,:),'linewidth',2)
    xlabel('Time','interpreter','latex')
    ylabel('$\beta_{eq}(t)$','interpreter','latex')
end

%% Get c(t) by solving the ODE from stochastic averaging.
disp("Solving the ODE to find c(t):")

ic = 0.00000001;
c = zeros(ndof, numel(time));
for i=1:ndof
    beta_eq_dof = beta_eq(i,:);
    omega_eq_2_dof = omega_eq_2(i,:);

    beta_eq_dof(1) = beta_eq_dof(2);
    omega_eq_2_dof(1) = omega_eq_2_dof(2);

    [t, c(i,:)] = ode45(@(t, c_aux) solve_c_mdof(t, c_aux, beta_eq_dof, omega_eq_2_dof, time), time, ic);
end

%% =====

for i=1:ndof

    for j=1:numel(time)
        t=time(j);
        sig2t = c(i,j);
        Sw = @(x)( evolutionary_power_spectrum(x, t) );
        Sx = @(x,y)( Sw(x)./( (omega_eq_2(i,j) - x.^2).^2 + (y*x).^2 )  );
        sfun = @(y) ( (sig2t - 2*integral(@(x)Sx(x,y),0,Inf)).^2  );

        bt = fminbnd(sfun,0.01,500);
        beq(i,j)=bt;

    end
    beq(i,1) = beq(i,2);

end

beta_eq = beq;

ic = 0.00000001;
c = zeros(ndof, numel(time));
for i=1:ndof
    beta_eq_dof = beta_eq(i,:);
    omega_eq_2_dof = omega_eq_2(i,:);

    beta_eq_dof(1) = beta_eq_dof(2);
    omega_eq_2_dof(1) = omega_eq_2_dof(2);

    [t, c(i,:)] = ode45(@(t, c_aux) solve_c_mdof(t, c_aux, beta_eq_dof, omega_eq_2_dof, time), time, ic);
end


%%
figure('color',[1 1 1]);
for i=1:ndof
    subplot(ndof,1,i); 
    hold on
    plot(time, varx_sl(i,:),'k-','linewidth',2) % SL
    plot(time, c(i,:)','r--','linewidth',2) % ODE
    if run_mcs
        plot(time_out,varx_mcs(i,:),'b:','linewidth',2) % MCS
    end
    legend('SL', 'SA', 'MCS','interpreter','latex')
    xlabel('Time','interpreter','latex')
    ylabel('$Var[x(t)]$','interpreter','latex')
end

%% Survival Probability

run_fps = input('Find the Survival Probability (Yes=1, No=0):');

if run_fps

    P=survival_probability(barrier,c,time,30,beta_eq,5);

    
    figure('color',[1 1 1]);
    for i=1:ndof
        fpt = first_passage_time(:,i);
        fpt = fpt(fpt>0);
        [fpp,tfp]=ksdensity(fpt,'width',0.001,'Function','survivor');
        subplot(ndof,1,i); 
        hold on
        scatter(fpt,linspace(0,1,numel(fpt)),10,'r','filled','MarkerFaceAlpha',0.2);
        plot(time, P(i,:)','k','linewidth',2);
        plot(tfp, fpp,'b--','linewidth',2);
        
        legend('MCS','Analytical','KDE','interpreter','latex')
        title('First passage time')
        xlabel('Time')
        ylabel('Survivor Propability')
        xlim([0 T])
        ylim([0 1])
    end

end
%%


