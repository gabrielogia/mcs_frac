clc
clear
close all

files = dir('data/');
epx = 0.20;
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
aux = sprintf('Oscillator equivalent natural frequency. Nonlinearity: %.2f, q = %.2f', epx, 0.5);
title(aux, 'FontSize', 14)
legend('DOF 1', 'DOF 2', 'DOF 3')
grid(1)

subplot(2,1,2)
plot(time, beta_eq', 'linewidth',2)
xlabel('Time (s)', 'interpreter','latex', 'FontSize', 14)
ylabel('$\beta_{eq}(t) $ (Ns/m)','interpreter','latex', 'FontSize', 14)
aux = sprintf('Oscillator equivalent damping. Nonlinearity: %.2f, q = %.2f', epx, 0.5);
title(aux, 'FontSize', 14)
legend('DOF 1', 'DOF 2', 'DOF 3')
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
            aux = sprintf('Oscillator equivalent natural frequency. DOF: %d; Nonlinearity: %.2f', dof, epx);
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
            aux = sprintf('Oscillator equivalent damping. DOF: %d; Nonlinearity: %.2f', dof, epx);
            title(aux, 'FontSize', 14)
            grid;
            legend('q = 0.35', 'q = 0.5', 'q = 0.75')
        end
    end
end

%% Plot first-passage PDF
str1 = 'pdfs_';

for i = 1:1:numel(files)
    if (contains(files(i).name, 'displacement_variance'))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(8));

        files(i).name

        c = vec.c;

        for i=1:ndof
            smaxt(i) = max(c(i,:));
        end
        
        smaxi = max(smaxt);
        for i=1:ndof
            barrier(i) = lam*sqrt(smaxi);
        end
    end
    
    if(contains(files(i).name, str1) && contains(files(i).name, str2))
        files(i).name
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(7));

        q = 0.5;
        if (str2double(str_split(7)) == q)
            av = vec.av;
            pa = vec.pa;
            pr = vec.pr;
            time_out = vec.time_out;

            bar = ones(size(time_out))*barrier(1);
            ha = ones(size(time_out))*1000;
    
            figure(4)
            colormap jet
            subplot(3,2,1)
            hold on
            surf(time_out,av,pr(:,:,3));
            plot3(time_out,bar,ha,'r','linewidth',2)
            clim([0,20])
            shading interp
            xlabel('Time (s)','interpreter','latex', 'FontSize', 14)
            ylabel('Amplitude (m)','interpreter','latex', 'FontSize', 14)
            aux = sprintf('Empirical probability density function. Fractional: %.2f; DOF: 3', vec(:).q);
            title(aux, 'FontSize', 14)
    
            subplot(3,2,2)
            hold on
            surf(time_out,av,pa(:,:,3));
            plot3(time_out,bar,ha,'r','linewidth',2)
            clim([0,20])
            shading interp
            xlabel('Time (s)','interpreter','latex', 'FontSize', 14)
            ylabel('Amplitude (m)','interpreter','latex', 'FontSize', 14)
            aux = sprintf('Analytical probability density function. Fractional: %.2f; DOF: 3', vec(:).q);
            title(aux, 'FontSize', 14)

            subplot(3,2,3)
            hold on
            surf(time_out,av,pr(:,:,2));
            clim([0,20])
            shading interp
            xlabel('Time (s)','interpreter','latex', 'FontSize', 14)
            ylabel('Amplitude (m)','interpreter','latex', 'FontSize', 14)
            aux = sprintf('Empirical probability density function. Fractional: %.2f; DOF: 2', vec(:).q);
            title(aux, 'FontSize', 14)
    
            subplot(3,2,4)
            hold on
            surf(time_out,av,pa(:,:,2));
            clim([0,20])
            shading interp
            xlabel('Time (s)','interpreter','latex', 'FontSize', 14)
            ylabel('Amplitude (m)','interpreter','latex', 'FontSize', 14)
            aux = sprintf('Analytical probability density function. Fractional: %.2f; DOF: 2', vec(:).q);
            title(aux, 'FontSize', 14)

            subplot(3,2,5)
            hold on
            surf(time_out,av,pr(:,:,1));
            clim([0,20])
            shading interp
            xlabel('Time (s)','interpreter','latex', 'FontSize', 14)
            ylabel('Amplitude (m)','interpreter','latex', 'FontSize', 14)
            aux = sprintf('Empirical probability density function. Fractional: %.2f; DOF: 1', vec(:).q);
            title(aux, 'FontSize', 14)
    
            subplot(3,2,6)
            hold on
            surf(time_out,av,pa(:,:,1));
            clim([0,20])
            shading interp
            xlabel('Time (s)','interpreter','latex', 'FontSize', 14)
            ylabel('Amplitude (m)','interpreter','latex', 'FontSize', 14)
            aux = sprintf('Analytical probability density function. Fractional: %.2f; DOF: 1', vec(:).q);
            title(aux, 'FontSize', 14)
            break
        end
    end
end

%% survival probability