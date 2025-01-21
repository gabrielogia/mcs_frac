clc
clear
close all

files = dir('data/');

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

set(groot,'defaultAxesFontSize',12)

%% Equivalent damping and natural frequencies for the same q
vec_omega = load('data/omegaeq_oscillator_duffing_ndof_3_fractional_0.50_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.25_powerspectrum_eps_S0_0.20_duffingparameter_epx_1.00');
vec_beta = load('data/betaeq_oscillator_duffing_ndof_3_fractional_0.50_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.25_powerspectrum_eps_S0_0.20_duffingparameter_epx_1.00');

time = vec_omega.time;
omega_eq_2 = vec_omega.omega_eq_2;
beta_eq = vec_beta.beta_eq;

fig = figure(1);
subplot(2,1,1)
plot(time, omega_eq_2(1,:)', '--', 'linewidth',2)
hold on
plot(time, omega_eq_2(2,:)', '-.', 'linewidth',2)
plot(time, omega_eq_2(3,:)', 'linewidth',2)
xlabel('Time (s)')
ylabel('$\omega_{eq}^2$  (N/m)')
xlim([0 11])
legend('DOF 1', 'DOF 2', 'DOF 3')
grid

subplot(2,1,2)
plot(time, beta_eq(1,:)', '--', 'linewidth',2)
hold on
plot(time, beta_eq(2,:)', '-.', 'linewidth',2)
plot(time, beta_eq(3,:)', 'linewidth',2)
xlabel('Time (s)')
ylabel('$\beta_{eq}$ (Ns/m)')
xlim([0 11])
legend('DOF 1', 'DOF 2', 'DOF 3', 'Location', 'southeast')
grid

set(fig,'papersize',[7 5.5]);
print(fig,'plots/equivalent_stiff_damp_same_q_same_e','-dpdf')

%% Equivalent stiffness and damping
epx = 1;
str1 = 'omegaeq_';
str2 = sprintf('epx_%.2f', epx);

markers = ["--", "-.", "-"];
letters = ["c", "b", "a"];

for i = 1:1:numel(files)
    if(contains(files(i).name, str1) && contains(files(i).name, str2))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(7));

        if (vec(:).q == 0.35)
            marker = markers(1);
        elseif (vec(:).q == 0.5)
            marker = markers(2);
        else
            marker = markers(3);
        end

        ndof = size(vec.omega_eq_2);
        ndof = ndof(1);

        time = vec.time;
        omega_eq_2 = vec.omega_eq_2;
    
        fig = figure(2);
        for dof = 1:1:ndof
            subplot(ndof,2, -2*dof + (2*ndof + 1)); 
            hold on
            plot(time,omega_eq_2(dof,:), marker, 'linewidth',2) % MCS
            xlabel('Time (s)')
            ylabel('$\omega^2_{eq} $ (N/m)')
            aux = sprintf("%s) DOF %d", letters(dof), dof);
            title(aux)
            grid;
            legend('q = 0.35', 'q = 0.5', 'q = 0.75')
            xlim([0 11])
        end
    end
end

str1 = 'betaeq_';

letters = ["f", "e", "d"];

for i = 1:1:numel(files)
    if(contains(files(i).name, str1) && contains(files(i).name, str2))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(7));

        if (vec(:).q == 0.35)
            marker = markers(1);
        elseif (vec(:).q == 0.5)
            marker = markers(2);
        else
            marker = markers(3);
        end

        ndof = size(vec.beta_eq);
        ndof = ndof(1);

        time = vec.time;
        beta_eq = vec.beta_eq;
    
        figure(2);
        for dof = 1:1:ndof
            subplot(ndof,2, -2*(dof - ndof - 1)); 
            hold on
            plot(time,beta_eq(dof,:), marker, 'linewidth',2) % MCS
            xlabel('Time (s)')
            ylabel('$\beta_{eq} $ (Ns/m)')
            aux = sprintf("%s) DOF %d", letters(dof), dof);
            title(aux)
            grid;
            if (dof == 3)
                legend('q = 0.35', 'q = 0.5', 'q = 0.75', 'Location', 'northwest')
            else
                legend('q = 0.35', 'q = 0.5', 'q = 0.75')
            end
            xlim([0 11])
        end
    end
end

set(fig,'papersize',[8.5 7.5], 'Position',[200 200 700 550]);
print(fig,'plots/equivalent_stiffness_and_beta_different_q','-dpdf')

%% Plot amplitude PDF
lam = 0.25;
str1 = 'pdfs_';
epx = 1;
ndof = 3;
str2 = sprintf('epx_%.2f', epx);
q = 0.50;

letters = ["c" "b" "a" "f" "e" "d"];

for i = 1:1:numel(files)
    if (contains(files(i).name, 'displacement_variance') && contains(files(i).name, str2))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(8));

        if (str2double(str_split(8)) == q)
            c = vec.c;
    
            for j=1:ndof
                smaxt(j) = max(c(j,:));
            end
            
            smaxi = max(smaxt);
            for j=1:ndof
                barrier(j) = lam*sqrt(smaxi);
            end
        end
    end
    
    if(contains(files(i).name, str1) && contains(files(i).name, str2))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(7));

        if (str2double(str_split(7)) == q)
            av = vec.av;
            pa = vec.pa;
            pr = vec.pr;
            time_out = vec.time_out;

            ha = ones(size(time_out))*1000;

            for j = 1:1:ndof
                fig = figure(4);
                colormap jet
                subplot(3,2,2*j-1)
                hold on
                surf(time_out,av,pr(:,:,3-j+1));
                plot3(time_out,ones(size(time_out))*barrier(3-j+1),ha,'r','linewidth',2)
                clim([0,150])
                ylim([0,max(av)])
                xlim([0 11])
                shading interp
                xlabel('Time (s)')
                ylabel('Amplitude (m)')
                aux = sprintf('%s) MCS PDF. DOF %d', letters(3-j+1), 3-j+1);
                title(aux)
        
                subplot(3,2,2*j)
                hold on
                surf(time_out,av,pa(:,:,3-j+1));
                plot3(time_out,ones(size(time_out))*barrier(3-j+1),ha,'r','linewidth',2)
                clim([0,150])
                ylim([0,max(av)])
                xlim([0 11])
                shading interp
                xlabel('Time (s)')
                ylabel('Amplitude (m)')
                aux = sprintf('%s) Approximate analytical PDF. DOF: %d', letters(3-j+4), 3-j+1);
                title(aux)
            end
            break
        end
    end
end

set(fig,'papersize',[9.5 7.5], 'Position',[200 200 800 550]);
print(fig,'plots/amplitude_pdf','-dpdf')

%% survival probability
str1 = 'firsttimepassage_';

k = 1;
for i = 1:1:numel(files)
    if (contains(files(i).name, str1) && contains(files(i).name, str2))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(7));

        P = vec.P;
        time = vec.time;
        fpp = vec.fpp;
        tfp = vec.tfp;
        T = 2;

        fig = figure(5);
        subplot(ndof,3,4-k); 
        hold on
        plot(time, P(3,:)','k','linewidth',2);
        plot(tfp, fpp,'r--','linewidth',2);
        
        legend('Analytical','MCS')
        aux = sprintf('$q = %.2f$; DOF: %d', vec(:).q, 3);
        title(aux)
        xlabel('Time (s)')
        ylabel('Survival propability')
        xlim([0 T])
        ylim([0 1])

        subplot(ndof,3,7-k); 
        hold on
        plot(time, P(2,:)','k','linewidth',2);
        plot(tfp, fpp,'r--','linewidth',2);
        
        legend('Analytical','MCS')
        aux = sprintf('$q = %.2f$; DOF: %d', vec(:).q, 2);
        title(aux)
        xlabel('Time (s)')
        ylabel('Survival propability')
        xlim([0 T])
        ylim([0 1])

        subplot(ndof,3,10-k); 
        hold on
        plot(time, P(1,:)','k','linewidth',2);
        plot(tfp, fpp,'r--','linewidth',2);
        
        legend('Analytical','MCS')
        aux = sprintf('$q = %.2f$; DOF: %d', vec(:).q, 1);
        title(aux)
        xlabel('Time (s)')
        ylabel('Survival propability')
        xlim([0 T])
        ylim([0 1])

        k = k + 1;
    end
end

set(fig,'papersize',[9.5 7.5], 'Position',[200 200 800 550]);
print(fig,'plots/survival_prop','-dpdf')

%% Equivalent stiffness and damping for diferent non-linearities
str0 = 'duffing';
str1 = 'omegaeq_';
q = 0.50;
str2 = sprintf('fractional_%.2f', q);
str3 = sprintf('mcssamples_%d', 14000);

markers = ["--", "-.", "-"];
letters = ["c", "b", "a"];

for i = 1:1:numel(files)
    if(contains(files(i).name, str0) && contains(files(i).name, str1) ...
            && contains(files(i).name, str2) && contains(files(i).name, str3))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        str_split = strsplit(string(str_split(24)), ".mat");
        vec(:).epx = str2double(str_split(1));

        if (vec(:).epx == 0.50)
            marker = markers(1);
        elseif (vec(:).epx == 1.00)
            marker = markers(2);
        else
            marker = markers(3);
        end

        ndof = size(vec.omega_eq_2);
        ndof = ndof(1);

        time = vec.time;
        omega_eq_2 = vec.omega_eq_2;
    
        fig = figure(4);
        for dof = 1:1:ndof
            subplot(ndof,2, -2*dof + (2*ndof + 1)); 
            hold on
            plot(time,omega_eq_2(dof,:), marker, 'linewidth',2) % MCS
            xlabel('Time (s)')
            ylabel('$\omega^2_{eq} $ (N/m)')
            aux = sprintf("%s) DOF %d", letters(dof), dof);
            title(aux)
            grid
            legend('$\epsilon$ = 0.50', '$\epsilon$ = 1.00', '$\epsilon$ = 2.00')
        end
    end
end

str1 = 'betaeq_';
letters = ["f", "e", "d"];

for i = 1:1:numel(files)
    if(contains(files(i).name, str0) && contains(files(i).name, str1) ...
            && contains(files(i).name, str2) && contains(files(i).name, str3))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        str_split = strsplit(string(str_split(24)), ".mat");
        vec(:).epx = str2double(str_split(1));

        if (vec(:).epx == 0.50)
            marker = markers(1);
        elseif (vec(:).epx == 1.00)
            marker = markers(2);
        else
            marker = markers(3);
        end

        ndof = size(vec.beta_eq);
        ndof = ndof(1);

        time = vec.time;
        beta_eq = vec.beta_eq;
    
        figure(4);
        for dof = 1:1:ndof
            subplot(ndof,2, -2*(dof - ndof - 1)); 
            hold on
            plot(time,beta_eq(dof,:), marker, 'linewidth',2) % MCS
            grid
            xlabel('Time (s)')
            ylabel('$\beta_{eq} $ (Ns/m)')
            aux = sprintf("%s) DOF %d", letters(dof), dof);
            title(aux)
            legend('$\epsilon$ = 0.50', '$\epsilon$ = 1.00', '$\epsilon$ = 2.00', 'Location', 'southeast')
        end
    end
end

set(fig,'papersize',[8.5 7.5], 'Position',[200 200 700 550]);
print(fig,'plots/equivalent_stiffines_and_damping_different_epi','-dpdf')

%% displacement same q and same ep
load('data/displacement_variance_oscillator_duffing_ndof_3_fractional_0.75_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.25_powerspectrum_eps_S0_0.20_duffingparameter_epx_1.00.mat')

fig = figure;
for i=1:ndof
    subplot(ndof,1,i); 
    hold on
    plot(time, sqrt(varx_sl(i,:)),'k-','linewidth',2) % SL
    plot(time, sqrt(c(i,:))','r--','linewidth',2) % ODE
    plot(time_out,sqrt(varx_mcs(i,:)),'b:','linewidth',2) % MCS
    legend('SL', 'SA', 'MCS', 'Location', 'southeast')
    xlabel('Time (s)')
    ylabel('$\sigma[x(t)] (m)$')
    xlim([0 11])
    aux = sprintf("DOF %d", i);
    title(aux)
end

set(fig,'papersize',[8.5 7.5], 'Position',[200 200 700 550]);
print(fig,'plots/displacement_std_same_q_same_epx','-dpdf')