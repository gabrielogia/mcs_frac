clc
clear
close all

files = dir('data/');
epx = 0.20;

%% Equivalent damping and natural frequencies for the same q
vec_omega = load(strcat('data/', files(i).name));
vec_beta = load('betaeq_oscillator_duffing_ndof_3_fractional_0.50_nonlinearity_0.90_dt_0.0010_mcssamples_12000_damping_20.00_stiffness_200.00');

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
    
        figure(1)
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
    
        figure(2)
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

%% Displacement variance
str1 = 'displacement_';

k = 1;
for i = 1:1:numel(files)
    if(contains(files(i).name, str1) && contains(files(i).name, str2))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(8));

        ndof = size(vec.varx_sl);
        ndof = ndof(1);

        time = vec.time;
        varx_sl = vec.varx_sl;
        c = vec.c;
        time_out = vec.time_out;
        varx_mcs = vec.varx_mcs;

        figure(3)
        for dof=1:1:ndof
            subplot(ndof,ndof,k); 
            hold on
            plot(time, varx_sl(dof,:),'k-','linewidth',2) % SL
            plot(time, c(dof,:)','r--','linewidth',2) % ODE
            plot(time_out,varx_mcs(dof,:),'b:','linewidth',2) % MCS
            legend('SL', 'SA', 'MCS','interpreter','latex')
            xlabel('Time (s)','interpreter','latex', 'FontSize', 14)
            ylabel('$Var[x(t)]$ (m)','interpreter','latex', 'FontSize', 14)
            aux = sprintf('Oscillator displacement variance. DOF: %d; Fractional: %.2f', dof, vec(:).q);
            title(aux, 'FontSize', 14)
            grid(1)
            k = k + 1;
        end
    end
end

%% Plot first-passage PDF
str1 = 'pdfs_';

for i = 1:1:numel(files)
    if(contains(files(i).name, str1) && contains(files(i).name, str2))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(7));

        q = 0.5;
        if (str2double(str_split(7)) == q)
            av = vec.av;
            pa = vec.pa;
            pr = vec.pr;
            time_out = vec.time_out;
    
            figure(4)
            colormap jet
            subplot(2,1,1)
            hold on
            surf(time_out,av,pr(:,:,2));
            clim([0,20])
            shading interp
            xlabel('Time (s)','interpreter','latex', 'FontSize', 14)
            ylabel('Amplitude (m)','interpreter','latex', 'FontSize', 14)
            aux = sprintf('Empirical probability density function. Fractional: %.2f', vec(:).q);
            title(aux, 'FontSize', 14)
    
            subplot(2,1,2)
            hold on
            surf(time_out,av,pa(:,:,2));
            clim([0,20])
            shading interp
            xlabel('Time (s)','interpreter','latex', 'FontSize', 14)
            ylabel('Amplitude (m)','interpreter','latex', 'FontSize', 14)
            aux = sprintf('Analytical probability density function. Fractional: %.2f', vec(:).q);
            title(aux, 'FontSize', 14)
            break
        end
    end
end

%% survival probability