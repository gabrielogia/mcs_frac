%% main.m
clc
clear
close all

%% Structure data
tic 

% For reproducibility:
rng(1111);

% Oscillator ('bw', 'duffing')
oscillator = "bw";

% Is Base motion / non-stationary (excitation):
is_base = false;
nonstat = true;

% Number of DOFs:
ndof = 1;

% Fractional derivative:
q = 0.5; 

% Nonlinearity parameter:
epx = 1.9*ones(1,ndof);

% Mass, damping, and stiffness vectors: 
mass = 1*ones(1,ndof); 
damping = 2*ones(1,ndof);
stiffness = 20*ones(1,ndof);

% Bouc-Wen parameters
a_bw = 0.2*ones(1, ndof);
A_bw = 1;
beta_bw = 0.5;
gamma_bw = 0.5;
n_bw = 1;
y0_bw = 0.1;

% Yielding displacement.
xy=y0_bw;

% Maximum time:
T = 12;

% Barrier:
lam = 0.7;

% Time increment for the Monte Carlo simulation.
dT = 0.0001; %dT = 0.0001;

% Construct matrices M, C, and K:
if (oscillator == "duffing")
    [M, C, K] = get_mck(mass, damping, stiffness, ndof);
elseif (oscillator == "bw")
    [M, C, K] = get_mck_bw(mass, damping, stiffness, a_bw, ndof, y0_bw);
end

% Maximum frequency of the power spectrum:
fmax_ps = 150; 

% Number of samples in the MCS:
ns = 1;

% Discretization in time and frequency for the Statistical Linearization:
ntime = 200;
nfreq = 1000;

% Run MCS:
run_mcs = true;

% Find the Survival Probability (Yes=1, No=0):
run_fps = true;

barrier = 0.0001*ones(ndof,1);

%% Monte Carlo Simulation
if run_mcs
    disp('Running MCS:')
    
    if (oscillator == "duffing")
    [varx_mcs, time_out, first_passage_time,response,velocity,amplitude] = ...
        monte_carlo(ns,M,C,K,epx,q,mass,damping,stiffness,fmax_ps,...
        nonstat, is_base,T,dT, barrier);
    elseif (oscillator == "bw")
        [varx_mcs, time_out, first_passage_time,response,velocity, z, amplitude] = ...
        monte_carlo_bw_new(ns,M,C,K,q,fmax_ps,nonstat,is_base, T, dT, barrier, ndof, A_bw, gamma_bw, beta_bw, xy);
    end
end

%% PLOT HYSTERESIS LOOP

Ka = K(1:ndof,1:ndof);
K1a = K(1:ndof,ndof+1:2*ndof);

for i=1:numel(time_out)

    z0 = z(:,i);
    x0 = response(:,i);

    S(:,i) = Ka*x0 + K1a*z0;


end

figure
hold on
for i=1:ndof
    subplot(1,ndof,i)
    plot(response(i,:),S(i,:))
end