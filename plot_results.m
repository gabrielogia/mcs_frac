clc
clear
close all

files = dir('data/');
lam = 0.7;

%% Equivalent damping and natural frequencies for the same q
vec_omega = load('data/omegaeq_oscillator_duffing_ndof_3_fractional_0.50_nonlinearity_0.20_dt_0.0010_mcssamples_12000_damping_20.00_stiffness_200.00');
vec_beta = load('data/betaeq_oscillator_duffing_ndof_3_fractional_0.50_nonlinearity_0.20_dt_0.0010_mcssamples_12000_damping_20.00_stiffness_200.00');

time = vec_omega.time;
omega_eq_2 = vec_omega.omega_eq_2;
beta_eq = vec_beta.beta_eq;

fig = figure(1);
subplot(2,1,1)
plot(time, omega_eq_2', 'linewidth',2)
xlabel('Time (s)', 'interpreter','latex', 'FontSize', 12)
ylabel('$\omega_{eq}^2$  (N/m)','interpreter','latex', 'FontSize', 12)
%aux = sprintf("System's equivalent stiffness coefficient");
%title(aux, 'FontSize', 14)
legend('DOF 1', 'DOF 2', 'DOF 3')
grid

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
set(gca,'TickLabelInterpreter','latex')

subplot(2,1,2)
plot(time, beta_eq', 'linewidth',2)
xlabel('Time (s)', 'interpreter','latex', 'FontSize', 12)
ylabel('$\beta_{eq}$ (Ns/m)','interpreter','latex', 'FontSize', 12)
%aux = sprintf("Oscillator' equivalent damping coefficient");
%title(aux, 'FontSize', 14)
legend('DOF 1', 'DOF 2', 'DOF 3', 'Location', 'southeast')
grid

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
set(gca,'TickLabelInterpreter','latex')

set(fig,'papersize',[6.75 5.0]);
print(fig,'plots/equivalent_stiff_damp_same_q_same_e','-dpdf')

%% Equivalent natural frequency
epx = 0.2;
str1 = 'omegaeq_';
str2 = sprintf('nonlinearity_%.2f', epx);

for i = 1:1:numel(files)
    if(contains(files(i).name, str1) && contains(files(i).name, str2))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(7));

        ndof = size(vec.omega_eq_2);
        ndof = ndof(1);

        time = vec.time;
        omega_eq_2 = vec.omega_eq_2;
    
        fig = figure(2);
        for dof = 1:1:ndof
            subplot(ndof,1,ndof - dof + 1); 
            hold on
            plot(time,omega_eq_2(dof,:),'linewidth',2) % MCS
            xlabel('Time (s)', 'interpreter','latex', 'FontSize', 12)
            ylabel('$\omega^2_{eq} $ (N/m)','interpreter','latex', 'FontSize', 12)
            aux = sprintf("DOF: %d", dof);
            title(aux, 'FontSize', 10)
            grid;
            legend('q = 0.35', 'q = 0.5', 'q = 0.75')
        end
    end
end

set(fig,'papersize',[6.75 5.0]);
print(fig,'plots/equivalent_stiff_different_q','-dpdf')

%% Equivalent damping
str1 = 'betaeq_';
str2 = sprintf('nonlinearity_%.2f', epx);

for i = 1:1:numel(files)
    if(contains(files(i).name, str1) && contains(files(i).name, str2))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(7));

        ndof = size(vec.beta_eq);
        ndof = ndof(1);

        time = vec.time;
        beta_eq = vec.beta_eq;
    
        fig = figure(3);
        for dof = 1:1:ndof
            subplot(ndof,1,ndof - dof + 1); 
            hold on
            plot(time,beta_eq(dof,:),'linewidth',2) % MCS
            xlabel('Time (s)', 'interpreter','latex', 'FontSize', 12)
            ylabel('$\beta_{eq} $ (Ns/m)','interpreter','latex', 'FontSize', 12)
            aux = sprintf("DOF: %d", dof);
            title(aux, 'FontSize', 10)
            grid;
            legend('q = 0.35', 'q = 0.5', 'q = 0.75')
        end
    end
end

set(fig,'papersize',[6.75 5.4]);
print(fig,'plots/equivalent_beta_different_q','-dpdf')

%% Plot amplitude PDF
str1 = 'pdfs_';
epx = 0.2;
ndof = 3;
str2 = sprintf('nonlinearity_%.2f', epx);
q = 0.50;

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
                xlabel('Time (s)','interpreter','latex', 'FontSize', 12)
                ylabel('Amplitude (m)','interpreter','latex', 'FontSize', 12)
                aux = sprintf('MCS PDF. DOF: %d', 3-j+1);
                title(aux, 'FontSize', 10)
        
                subplot(3,2,2*j)
                hold on
                surf(time_out,av,pa(:,:,3-j+1));
                plot3(time_out,ones(size(time_out))*barrier(3-j+1),ha,'r','linewidth',2)
                clim([0,25])
                ylim([0,0.4])
                shading interp
                xlabel('Time (s)','interpreter','latex', 'FontSize', 12)
                ylabel('Amplitude (m)','interpreter','latex', 'FontSize', 12)
                aux = sprintf('Approximate analytical PDF. DOF: %d', 3-j+1);
                title(aux, 'FontSize', 10)
            end
            break
        end
    end
end

set(fig,'papersize',[6.75 5.4]);
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
        T = 20;

        fig = figure(5);
        subplot(ndof,3,k); 
        hold on
        plot(time, P(3,:)','k','linewidth',2);
        plot(tfp, fpp,'r--','linewidth',2);
        
        legend('Analytical','MCS','interpreter','latex')
        aux = sprintf('$q = %.2f$; DOF: %d', vec(:).q, 3);
        title(aux, 'FontSize', 13, 'interpreter','latex')
        xlabel('Time (s)','interpreter','latex', 'FontSize', 12)
        ylabel('Survival propability','interpreter','latex', 'FontSize', 12)
        xlim([0 T])
        ylim([0 1])

        subplot(ndof,3,k+3); 
        hold on
        plot(time, P(2,:)','k','linewidth',2);
        plot(tfp, fpp,'r--','linewidth',2);
        
        legend('Analytical','MCS','interpreter','latex')
        aux = sprintf('$q = %.2f$; DOF: %d', vec(:).q, 2);
        title(aux, 'FontSize', 13, 'interpreter','latex')
        xlabel('Time (s)','interpreter','latex', 'FontSize', 12)
        ylabel('Survival propability','interpreter','latex', 'FontSize', 12)
        xlim([0 T])
        ylim([0 1])

        subplot(ndof,3,k+6); 
        hold on
        plot(time, P(1,:)','k','linewidth',2);
        plot(tfp, fpp,'r--','linewidth',2);
        
        legend('Analytical','MCS','interpreter','latex')
        aux = sprintf('$q = %.2f$; DOF: %d', vec(:).q, 1);
        title(aux, 'FontSize', 13, 'interpreter','latex')
        xlabel('Time (s)','interpreter','latex', 'FontSize', 12)
        ylabel('Survival propability','interpreter','latex', 'FontSize', 12)
        xlim([0 T])
        ylim([0 1])

        k = k + 1;
    end
end

set(fig,'papersize',[6.8 5.5]);
print(fig,'plots/survival_prop','-dpdf')

%% Equivalent natural frequency for diferent non-linearities
str1 = 'omegaeq_';
q = 0.50;
str2 = sprintf('fractional_%.2f', q);

for i = 1:1:numel(files)
    if(contains(files(i).name, str1) && contains(files(i).name, str2))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(7));

        ndof = size(vec.omega_eq_2);
        ndof = ndof(1);

        time = vec.time;
        omega_eq_2 = vec.omega_eq_2;
    
        fig = figure(6);
        for dof = 1:1:ndof
            subplot(ndof,1,ndof - dof + 1); 
            hold on
            plot(time,omega_eq_2(dof,:),'linewidth',2)
            grid(1);
            xlabel('Time (s)', 'interpreter','latex', 'FontSize', 12)
            ylabel('$\omega^2_{eq}$ (N/m)','interpreter','latex', 'FontSize', 12)
            aux = sprintf('DOF: %d', dof);
            title(aux, 'FontSize', 10)
            legend('\epsilon = 0.20', '\epsilon = 0.90', '\epsilon = 1.40', '\epsilon = 2.20')
        end
    end
end

set(fig,'papersize',[7 5.3]);
print(fig,'plots/equivalent_stiff_different_epi','-dpdf')

%% Equivalent damping for diferent non-linearities
str1 = 'betaeq_';
q = 0.50;
str2 = sprintf('fractional_%.2f', q);

for i = 1:1:numel(files)
    if(contains(files(i).name, str1) && contains(files(i).name, str2))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(7));

        ndof = size(vec.beta_eq);
        ndof = ndof(1);

        time = vec.time;
        beta_eq = vec.beta_eq;
    
        fig = figure(7);
        for dof = 1:1:ndof
            subplot(ndof,1,ndof - dof + 1); 
            hold on
            plot(time,beta_eq(dof,:),'linewidth',2) % MCS
            grid(1);
            xlabel('Time (s)', 'interpreter','latex', 'FontSize', 12)
            ylabel('$\beta_{eq}$ (Ns/m)','interpreter','latex', 'FontSize', 12)
            aux = sprintf('DOF: %d', dof);
            title(aux, 'FontSize', 10)
            legend('\epsilon = 0.20', '\epsilon = 0.90', '\epsilon = 1.40', '\epsilon = 2.20', 'Location', 'southeast')
        end
    end
end

set(fig,'papersize',[7 5.3]);
print(fig,'plots/equivalent_damping_different_epi','-dpdf')

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
