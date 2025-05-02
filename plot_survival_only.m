%% main.m
clc
clear
close all
warning('off','all')

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

set(groot,'defaultAxesFontSize', 11)

%%

% Power spectrum (wn, eps)
power_spectrum = "eps";

% Oscillator ('bw', 'duffing')
oscillator = "soft";

% Number of DOFs:
ndof = 3;

% Power spectrum mag
S0 = 0.2;

% Fractional derivative:
q = 0.75; 

% Nonlinearity parameter:
epx = 1.0*ones(1,ndof);
alpha = 0.5*ones(1,ndof);

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
y0_bw = 0.001;

% Yielding displacement.
xy=y0_bw;

% Barrier:
lam = 0.25;

% Time increment for the Monte Carlo simulation.
dT = 1e-3;

% Number of samples in the MCS:
ns = 14000;

%% plot survival for different q and lambda
col = ["a" "b" "c"];

Q = [0.5];
numQ = numel(Q);

for jj=1:numQ
    q = Q(jj);
    str = sprintf('oscillator_%s_ndof_%d_fractional_%.2f_dt_%.4f_mcssamples_%d_damping_%.2f_stiffness_%.2f_barrier_%.2f_powerspectrum_%s_S0_%.2f', ...
        oscillator, ndof, q, dT, ns, max(damping), max(stiffness), 0.25, power_spectrum, S0);
    
    if (oscillator == "bw")
        str = strcat(str, sprintf('_bwparameters_a_%.2f_A_%.2f_beta_%.2f_gamma_%.2f_xy_%.2f', max(a_bw), A_bw, beta_bw, gamma_bw, xy));
    elseif (oscillator == "duffing")
        str = strcat(str, sprintf('_duffingparameter_epx_%.2f', max(epx)));
    else
        str = strcat(str, sprintf('_softparameter_alpha_%.2f', max(alpha)));
    end
    
    % load mcs previous results
    if exist(strcat('data/mcs/mcs_', str, '.mat'))
        load(strcat('data/mcs/mcs_', str, '.mat'))
        run_mcs = false;
    else
        run_mcs = true;
    end
    
    load(strcat('data/firsttimepassage_', str, '.mat'))
    load(strcat('data/omegaeq_', str, '.mat'))
    load(strcat('data/betaeq_', str, '.mat'))
    
    lam = [0.25, 0.50, 0.75];
    numB = numel(lam);

    for ii=1:numB
        count = 1;
        B = sqrt(smaxi)*lam(ii);
        tf = zeros(ns,ndof);

        for i=1:ns
            for j=1:ndof
                amp=amplitude(j,:,i);
                time_aux = time_out(abs(amp) > B);
            
                if numel(time_aux)==0
                    tf(i,j) = NaN;
                else
                    tf(i,j) = time_aux(1);
                end
            end
        end

        P = survival_probability(B*ones(ndof), c, time, numel(time), beta_eq, omega_eq_2, 15, S0, oscillator, q);
    
        fig = figure(1);
        for i=1:ndof
            fpt = tf(:,i);
            fpt = fpt(fpt>0);
            [fpp,tfp]=ksdensity(fpt,'width',0.1,'Function','survivor');
            subplot(ndof,numQ,i + (jj-1)*ndof);
            hold on
            if (lam(ii) == 0.25)
                plot(time, P(i,:)','k','linewidth',2);
                plot(tfp, fpp,'--g','linewidth',2);
            elseif (lam(ii) == 0.50)
                plot(time, P(i,:)','m','linewidth',2);
                plot(tfp, fpp,'--c','linewidth',2);
            else
                plot(time, P(i,:)','r','linewidth',2);
                plot(tfp, fpp,'--b','linewidth',2);
            end
            fprintf("MSE: %.4f, q = %.2f, lambda = %.2f, dof: %d\n", mse(P(i,1:numel(fpp)),fpp), q, lam(ii), i)
            title(sprintf("%s) DOF: %d", col(i + (jj-1)*ndof), i))
            xlabel('Time')
            ylabel('Survival')
            xlim([0 2.5])
            ylim([0 1])
            xticks([0 0.5 1 1.5 2 2.5])
            grid on;
        end
    end
end
%%
legend('Analytical: $\lambda$ = 0.25', 'MCS: $\lambda$ = 0.25', ...
    'Analytical: $\lambda$ = 0.50', 'MCS: $\lambda$ = 0.50', ...
    'Analytical: $\lambda$ = 0.75', 'MCS: $\lambda$ = 0.75',  ...
    'Orientation', 'horizontal')
%set(fig,'papersize',[6.0 5.5], 'Position',[200 200 1400 900]);
%print(fig,strcat('plots/survival_prop_only_', oscillator),'-dpng','-r1000')
