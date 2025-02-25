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

set(groot,'defaultAxesFontSize', 16)

%%
% Power spectrum (wn, eps)
power_spectrum = "eps";

% Oscillator ('bw', 'duffing')
oscillator = "bw";

% Number of DOFs:
ndof = 3;

% Power spectrum mag
S0 = 0.2;

% Fractional derivative:
q = 0.75; 

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
y0_bw = 0.001;

% Yielding displacement.
xy=y0_bw;

% Barrier:
lam = 0.25;

%% plot survival for different q and lambda
col = ["a", "b", "c"];
Q = [0.5 0.75 0.999];
numQ = numel(Q);

for jj=1:numQ
    q = Q(jj);
    str = sprintf('oscillator_%s_ndof_%d_fractional_%.2f_dt_%.4f_mcssamples_%d_damping_%.2f_stiffness_%.2f_barrier_%.2f_powerspectrum_%s_S0_%.2f', ...
        oscillator, ndof, q, 1e-3, 14000, max(damping), max(stiffness), 0.25, power_spectrum, S0);
    
    if (oscillator == "bw")
        str = strcat(str, sprintf('_bwparameters_a_%.2f_A_%.2f_beta_%.2f_gamma_%.2f_xy_%.2f', max(a_bw), A_bw, beta_bw, gamma_bw, xy));
    else
        str = strcat(str, sprintf('_duffingparameter_epx_%.2f', max(epx)));
    end
    
    load(strcat('data/firsttimepassage_', str, '.mat'))
    load(strcat('data/omegaeq_', str, '.mat'))
    load(strcat('data/betaeq_', str, '.mat'))
    
    lam = [0.25, 0.50, 0.75];
    numB = numel(lam);

    for ii=1:numB
        count = 1;
        B = sqrt(smaxi)*lam(ii);
        P = survival_probability(B*ones(ndof), c, time, numel(time), beta_eq, omega_eq_2, 15, S0);
    
        fig = figure(1);
        for i=1:ndof
            subplot(1,ndof+1,i);
            hold on
            if (lam(ii) == 0.25)
                if (q == 0.75)
                    plot(time, P(i,:)','k','linewidth',2);
                elseif (q == 0.50)
                    plot(time, P(i,:)','--k','linewidth',2);
                else
                    plot(time, P(i,:)',':k','linewidth',2);
                end
            elseif (lam(ii) == 0.50)
                if (q == 0.75)
                    plot(time, P(i,:)','b','linewidth',2);
                elseif (q == 0.50)
                    plot(time, P(i,:)','--b','linewidth',2);
                else
                    plot(time, P(i,:)',':b','linewidth',2);
                end
            else
                if (q == 0.75)
                    plot(time, P(i,:)','r','linewidth',2);
                elseif (q == 0.50)
                    plot(time, P(i,:)','--r','linewidth',2);
                else
                    plot(time, P(i,:)',':r','linewidth',2);
                end
            end
            title(sprintf("%s) DOF: %d", col(i), i))
            xlabel('Time')
            ylabel('Survival propability')
            xlim([0 4])
            ylim([0 1])
            xticks([0 1.0 2.0 3.0 4.0])
            grid on;
        end
    end
end
%%
lg = legend("$\lambda$ = 0.25, $q$ = 0.50", ...
       "$\lambda$ = 0.50, $q$ = 0.50", ...
       "$\lambda$ = 0.75, $q$ = 0.50", ...
       "$\lambda$ = 0.25, $q$ = 0.75", ...
       "$\lambda$ = 0.50, $q$ = 0.75", ...
       "$\lambda$ = 0.75, $q$ = 0.75", ...
       "$\lambda$ = 0.25, $q$ = 1.00", ...
       "$\lambda$ = 0.50, $q$ = 1.00", ...
       "$\lambda$ = 0.75, $q$ = 1.00");
lg.Position(1:2) = [0.68 0.46];
set(fig,'papersize',[6.0 5.5], 'Position',[200 200 900 250]);
print(fig,strcat('plots/survival_prop_only_', oscillator),'-dpng','-r1000')
