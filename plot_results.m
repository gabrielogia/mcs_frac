clc
clear
close all

files = dir('data/');
lam = 0.7;

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

set(groot,'defaultAxesFontSize',12)

%% Equivalent damping and natural frequencies for the same q
vec_omega = load('data/omegaeq_oscillator_duffing_ndof_3_fractional_0.50_nonlinearity_0.20_dt_0.0010_mcssamples_12000_damping_20.00_stiffness_200.00');
vec_beta = load('data/betaeq_oscillator_duffing_ndof_3_fractional_0.50_nonlinearity_0.20_dt_0.0010_mcssamples_12000_damping_20.00_stiffness_200.00');

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
%aux = sprintf("System's equivalent stiffness coefficient");
%title(aux, 'FontSize', 14)
legend('DOF 1', 'DOF 2', 'DOF 3')
grid

subplot(2,1,2)
plot(time, beta_eq(1,:)', '--', 'linewidth',2)
hold on
plot(time, beta_eq(2,:)', '-.', 'linewidth',2)
plot(time, beta_eq(3,:)', 'linewidth',2)
xlabel('Time (s)')
ylabel('$\beta_{eq}$ (Ns/m)')
%aux = sprintf("Oscillator' equivalent damping coefficient");
%title(aux, 'FontSize', 14)
legend('DOF 1', 'DOF 2', 'DOF 3', 'Location', 'southeast')
grid

set(fig,'papersize',[7 5.5]);
print(fig,'plots/equivalent_stiff_damp_same_q_same_e','-dpdf')

%% Equivalent stiffness and damping
epx = 0.2;
str1 = 'omegaeq_';
str2 = sprintf('nonlinearity_%.2f', epx);

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
        end
    end
end

str1 = 'betaeq_';
str2 = sprintf('nonlinearity_%.2f', epx);

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
        end
    end
end

set(fig,'papersize',[8.5 7.5], 'Position',[200 200 700 550]);
print(fig,'plots/equivalent_stiffness_and_beta_different_q','-dpdf')

%% Plot amplitude PDF
str1 = 'pdfs_';
epx = 0.2;
ndof = 3;
str2 = sprintf('nonlinearity_%.2f', epx);
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
                clim([0,25])
                ylim([0,0.4])
                shading interp
                xlabel('Time (s)')
                ylabel('Amplitude (m)')
                aux = sprintf('%s) MCS PDF. DOF %d', letters(3-j+1), 3-j+1);
                title(aux)
        
                subplot(3,2,2*j)
                hold on
                surf(time_out,av,pa(:,:,3-j+1));
                plot3(time_out,ones(size(time_out))*barrier(3-j+1),ha,'r','linewidth',2)
                clim([0,25])
                ylim([0,0.4])
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
close all
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
        T = 5;

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
str3 = sprintf('mcssamples_%d', 12000);

markers = ["--", "-.", "-"];
letters = ["c", "b", "a"];

for i = 1:1:numel(files)
    if(contains(files(i).name, str0) && contains(files(i).name, str1) ...
            && contains(files(i).name, str2) && contains(files(i).name, str3))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).epx = str2double(str_split(9));

        if vec(:).epx == 0.2
            continue
        end

        if (vec(:).epx == 0.9)
            marker = markers(1);
        elseif (vec(:).epx == 1.4)
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
            legend('$\epsilon$ = 0.90', '$\epsilon$ = 1.40', '$\epsilon$ = 2.20')
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
        vec(:).epx = str2double(str_split(9));

        if round(vec(:).epx,1) == 0.2
            continue
        end

        if (vec(:).epx == 0.9)
            marker = markers(1);
        elseif (vec(:).epx == 1.4)
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
            legend('$\epsilon$ = 0.90', '$\epsilon$ = 1.40', '$\epsilon$ = 2.20', 'Location', 'southeast')
        end
    end
end

set(fig,'papersize',[8.5 7.5], 'Position',[200 200 700 550]);
print(fig,'plots/equivalent_stiffines_and_damping_different_epi','-dpdf')

%%
a = load('data/omegaeq_oscillator_duffing_ndof_3_fractional_0.50_nonlinearity_0.20_dt_0.0010_mcssamples_12000_damping_20.00_stiffness_200.00.mat');
b = load('data/omegaeq_oscillator_duffing_ndof_3_fractional_0.50_nonlinearity_0.90_dt_0.0010_mcssamples_12000_damping_20.00_stiffness_200.00.mat');
c = load('data/omegaeq_oscillator_duffing_ndof_3_fractional_0.50_nonlinearity_1.40_dt_0.0010_mcssamples_12000_damping_20.00_stiffness_200.00.mat');
q = 0.5;

figure;
subplot(3,1,1)
hold on
dof = 3;
plot(a.time, a.omega_eq_2(dof,:),'linewidth',2);
plot(a.time, b.omega_eq_2(dof,:),'linewidth',2);
plot(a.time, c.omega_eq_2(dof,:),'linewidth',2);
xlabel('Time (s)', 'interpreter','latex', 'FontSize', 14)
ylabel('$\omega^2_{eq}(t) $ (rad/s)','interpreter','latex', 'FontSize', 14)
legend('$\epsilon = 0.20$', '$\epsilon = 0.90$', '$\epsilon = 1.40$', 'interpreter','latex', 'FontSize', 13)
aux = sprintf('Oscillator equivalent natural frequency. DOF: %d; Fractional: %.2f', dof, q);
title(aux, 'FontSize', 14)
grid;

subplot(3,1,2)
hold on
dof = 2;
plot(a.time, a.omega_eq_2(dof,:),'linewidth',2);
plot(a.time, b.omega_eq_2(dof,:),'linewidth',2);
plot(a.time, c.omega_eq_2(dof,:),'linewidth',2);
xlabel('Time (s)', 'interpreter','latex', 'FontSize', 14)
ylabel('$\omega^2_{eq}(t) $ (rad/s)','interpreter','latex', 'FontSize', 14)
legend('$\epsilon = 0.20$', '$\epsilon = 0.90$', '$\epsilon = 1.40$', 'interpreter','latex', 'FontSize', 13)
aux = sprintf('Oscillator equivalent natural frequency. DOF: %d; Fractional: %.2f', dof, q);
title(aux, 'FontSize', 14)
grid;

subplot(3,1,3)
hold on
dof = 1;
plot(a.time, a.omega_eq_2(dof,:),'linewidth',2);
plot(a.time, b.omega_eq_2(dof,:),'linewidth',2);
plot(a.time, c.omega_eq_2(dof,:),'linewidth',2);
xlabel('Time (s)', 'interpreter','latex', 'FontSize', 14)
ylabel('$\omega^2_{eq}(t) $ (rad/s)','interpreter','latex', 'FontSize', 14)
legend('$\epsilon = 0.20$', '$\epsilon = 0.90$', '$\epsilon = 1.40$', 'interpreter','latex', 'FontSize', 13)
aux = sprintf('Oscillator equivalent natural frequency. DOF: %d; Fractional: %.2f', dof, q);
title(aux, 'FontSize', 14)
grid;
