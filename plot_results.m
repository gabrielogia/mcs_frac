clc
clear
close all

files = dir('data/');
epx = 0.2;
lam = 0.7;

%% Equivalent damping and natural frequencies for the same q
vec_omega = load('data/omegaeq_oscillator_duffing_ndof_3_fractional_0.50_nonlinearity_0.20_dt_0.0010_mcssamples_12000_damping_20.00_stiffness_200.00');
vec_beta = load('data/betaeq_oscillator_duffing_ndof_3_fractional_0.50_nonlinearity_0.20_dt_0.0010_mcssamples_12000_damping_20.00_stiffness_200.00');

time = vec_omega.time;
omega_eq_2 = vec_omega.omega_eq_2;
beta_eq = vec_beta.beta_eq;

figure(1)
subplot(2,1,1)
plot(time, omega_eq_2', 'linewidth',2)
xlabel('Time (s)', 'interpreter','latex', 'FontSize', 14)
ylabel('$\omega^2_{eq}(t) $ (rad/s)','interpreter','latex', 'FontSize', 14)
aux = sprintf("Oscillator' equivalent natural frequencies");
title(aux, 'FontSize', 14)
legend('DOF 1', 'DOF 2', 'DOF 3')
grid(1)

subplot(2,1,2)
plot(time, beta_eq', 'linewidth',2)
xlabel('Time (s)', 'interpreter','latex', 'FontSize', 14)
ylabel('$\beta_{eq}(t) $ (Ns/m)','interpreter','latex', 'FontSize', 14)
aux = sprintf("Oscillator' equivalent damping");
title(aux, 'FontSize', 14)
legend('DOF 1', 'DOF 2', 'DOF 3', 'Location', 'southeast')
grid(1)

%% Equivalent natural frequency
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
    
        figure(2)
        for dof = 1:1:ndof
            subplot(ndof,1,ndof - dof + 1); 
            hold on
            plot(time,omega_eq_2(dof,:),'linewidth',2) % MCS
            xlabel('Time (s)', 'interpreter','latex', 'FontSize', 14)
            ylabel('$\omega^2_{eq}(t) $ (rad/s)','interpreter','latex', 'FontSize', 14)
            aux = sprintf("Oscillator' equivalent natural frequency. DOF: %d;", dof);
            title(aux, 'FontSize', 14)
            grid;
            legend('q = 0.35', 'q = 0.5', 'q = 0.75')
        end
    end
end

%% Equivalent damping
str1 = 'betaeq_';

for i = 1:1:numel(files)
    if(contains(files(i).name, str1) && contains(files(i).name, str2))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(7));

        ndof = size(vec.beta_eq);
        ndof = ndof(1);

        time = vec.time;
        beta_eq = vec.beta_eq;
    
        figure(3)
        for dof = 1:1:ndof
            subplot(ndof,1,ndof - dof + 1); 
            hold on
            plot(time,beta_eq(dof,:),'linewidth',2) % MCS
            xlabel('Time (s)', 'interpreter','latex', 'FontSize', 14)
            ylabel('$\beta_{eq}(t) $ (Ns/m)','interpreter','latex', 'FontSize', 14)
            aux = sprintf("Oscillator equivalent damping. DOF: %d;", dof);
            title(aux, 'FontSize', 14)
            grid;
            legend('q = 0.35', 'q = 0.5', 'q = 0.75')
        end
    end
end

%% Plot amplitude PDF
str1 = 'pdfs_';
q = 0.75;

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
                figure(4)
                colormap jet
                subplot(3,2,2*j-1)
                colorbar
                hold on
                surf(time_out,av,pr(:,:,3-j+1));
                plot3(time_out,ones(size(time_out))*barrier(3-j+1),ha,'r','linewidth',2)
                clim([0,35])
                ylim([0,0.35])
                shading interp
                xlabel('Time (s)','interpreter','latex', 'FontSize', 14)
                ylabel('Amplitude (m)','interpreter','latex', 'FontSize', 14)
                aux = sprintf('MCS probability density function. DOF: %d', 3-j+1);
                title(aux, 'FontSize', 14)
        
                subplot(3,2,2*j)
                hold on
                colorbar
                surf(time_out,av,pa(:,:,3-j+1));
                plot3(time_out,ones(size(time_out))*barrier(3-j+1),ha,'r','linewidth',2)
                clim([0,35])
                ylim([0,0.35])
                shading interp
                xlabel('Time (s)','interpreter','latex', 'FontSize', 14)
                ylabel('Amplitude (m)','interpreter','latex', 'FontSize', 14)
                aux = sprintf('Analytical probability density function. DOF: %d', 3-j+1);
                title(aux, 'FontSize', 14)
            end
            break
        end
    end
end

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
        T = 20;

        figure(5)
        subplot(ndof,3,k); 
        hold on
        plot(time, P(3,:)','k','linewidth',2);
        plot(tfp, fpp,'r--','linewidth',2);
        
        legend('Analytical','MCS','interpreter','latex')
        aux = sprintf('Survival probability. Fractional: %.2f; DOF: %d', vec(:).q, 3);
        title(aux, 'FontSize', 14)
        xlabel('Time','interpreter','latex', 'FontSize', 14)
        ylabel('Propability','interpreter','latex', 'FontSize', 14)
        xlim([0 T])
        ylim([0 1])

        subplot(ndof,3,k+3); 
        hold on
        plot(time, P(2,:)','k','linewidth',2);
        plot(tfp, fpp,'r--','linewidth',2);
        
        legend('Analytical','MCS','interpreter','latex')
        aux = sprintf('Survival probability. Fractional: %.2f; DOF: %d', vec(:).q, 2);
        title(aux, 'FontSize', 14)
        xlabel('Time','interpreter','latex', 'FontSize', 14)
        ylabel('Propability','interpreter','latex', 'FontSize', 14)
        xlim([0 T])
        ylim([0 1])

        subplot(ndof,3,k+6); 
        hold on
        plot(time, P(1,:)','k','linewidth',2);
        plot(tfp, fpp,'r--','linewidth',2);
        
        legend('Analytical','MCS','interpreter','latex')
        aux = sprintf('Survival probability. Fractional: %.2f; DOF: %d', vec(:).q, 1);
        title(aux, 'FontSize', 14)
        xlabel('Time','interpreter','latex', 'FontSize', 14)
        ylabel('Propability','interpreter','latex', 'FontSize', 14)
        xlim([0 T])
        ylim([0 1])

        k = k + 1;
    end
end

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
    
        figure(6)
        for dof = 1:1:ndof
            subplot(ndof,1,ndof - dof + 1); 
            hold on
            plot(time,omega_eq_2(dof,:),'linewidth',2)
            xlabel('Time (s)', 'interpreter','latex', 'FontSize', 14)
            ylabel('$\omega^2_{eq}(t) $ (rad/s)','interpreter','latex', 'FontSize', 14)
            aux = sprintf('Oscillator equivalent natural frequency. DOF: %d; Fractional: %.2f', dof, q);
            title(aux, 'FontSize', 14)
            grid;
            legend('$\epsilon = 0.20$', '$\epsilon = 0.90$', '$\epsilon = 1.40$', 'interpreter','latex', 'FontSize', 12)
        end
    end
end

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
    
        figure(7)
        for dof = 1:1:ndof
            subplot(ndof,1,ndof - dof + 1); 
            hold on
            plot(time,beta_eq(dof,:),'linewidth',2) % MCS
            xlabel('Time (s)', 'interpreter','latex', 'FontSize', 14)
            ylabel('$\beta_{eq}(t) $ (Ns/m)','interpreter','latex', 'FontSize', 14)
            aux = sprintf('Oscillator equivalent damping. DOF: %d; Fractional: %.2f', dof, q);
            title(aux, 'FontSize', 14)
            grid;
            legend('$\epsilon = 0.20$', '$\epsilon = 0.90$', '$\epsilon = 1.40$', 'interpreter','latex', 'FontSize', 13,'Location','southeast')
        end
    end
end

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