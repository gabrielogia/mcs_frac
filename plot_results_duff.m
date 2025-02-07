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

set(groot,'defaultAxesFontSize',16)

%% Equivalent damping and natural frequencies for the same q
vec_omega = load('data/omegaeq_oscillator_duffing_ndof_3_fractional_0.75_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.25_powerspectrum_eps_S0_0.20_duffingparameter_epx_1.00');
vec_beta = load('data/betaeq_oscillator_duffing_ndof_3_fractional_0.75_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.25_powerspectrum_eps_S0_0.20_duffingparameter_epx_1.00');

time = vec_omega.time;
omega_eq_2 = vec_omega.omega_eq_2;
beta_eq = vec_beta.beta_eq;

fig = figure(1);
subplot(1,2,1)
plot(time, omega_eq_2(1,:)', '--', 'linewidth',2)
hold on
plot(time, omega_eq_2(2,:)', '-.', 'linewidth',2)
plot(time, omega_eq_2(3,:)', 'linewidth',2)
xlabel('Time (s)')
ylabel('$\omega_{e}^2$  (N/m)')
xlim([0 4])
title("a) Effective time-varying stiffness")
legend('DOF 1', 'DOF 2', 'DOF 3')
grid

subplot(1,2,2)
plot(time, beta_eq(1,:)', '--', 'linewidth',2)
hold on
plot(time, beta_eq(2,:)', '-.', 'linewidth',2)
plot(time, beta_eq(3,:)', 'linewidth',2)
xlabel('Time (s)')
ylabel('$\beta_{e}$ (Ns/m)')
xlim([0 4])
title("b) Effective time-varying damping")
legend('DOF 1', 'DOF 2', 'DOF 3')
grid

set(fig,'papersize',[6.0 5.5], 'Position', [200 200 900 350]);
print(fig,'plots/equivalent_stiff_damp_same_q_same_e','-dpng','-r1000')

%% Equivalent stiffness and damping
epx = 1;
str1 = 'omegaeq_';
str2 = sprintf('epx_%.2f', epx);

markers = ["-", "-."];
letters = ["a", "b", "b"];

for i = 1:1:numel(files)
    if(contains(files(i).name, str1) && contains(files(i).name, str2))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(7));

        if (vec(:).q == 0.50)
            marker = markers(1);
        elseif (vec(:).q == 0.75)
            marker = markers(2);
        else
            marker = "";
        end

        ndof = size(vec.omega_eq_2);
        ndof = ndof(1);

        time = vec.time;
        omega_eq_2 = vec.omega_eq_2;
    
        fig = figure(2);
        if (marker ~= "")
            for dof = 1:1:ndof
                subplot(2,ndof, dof); 
                hold on
                plot(time,omega_eq_2(dof,:), marker, 'linewidth',2)
                xlabel('Time (s)')
                ylabel('$\omega^2_{e} $ (N/m)')
                aux = sprintf("%s) DOF %d", letters(dof), dof);
                title(aux)
                grid(1);
                xlim([0 4])
                ylim([0 1500])
                legend('q = 0.5', 'q = 0.75')
            end
        end
    end
end

str1 = 'betaeq_';
letters = ["d", "e", "f"];

for i = 1:1:numel(files)
    if(contains(files(i).name, str1) && contains(files(i).name, str2))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(7));

        if (vec(:).q == 0.50)
            marker = markers(1);
        elseif (vec(:).q == 0.75)
            marker = markers(2);
        else
            marker = "";
        end

        ndof = size(vec.beta_eq);
        ndof = ndof(1);

        time = vec.time;
        beta_eq = vec.beta_eq;
    
        figure(2);
        if (marker ~= "")
            for dof = 1:1:ndof
                subplot(2,ndof, ndof+dof); 
                hold on
                plot(time,beta_eq(dof,:), marker, 'linewidth',2)
                xlabel('Time (s)')
                ylabel('$\beta_{e} $ (Ns/m)')
                aux = sprintf("%s) DOF %d", letters(dof), dof);
                title(aux)
                grid(1);
                xlim([0 4])
                ylim([0 40])
                legend('q = 0.5', 'q = 0.75')
            end
        end
    end
end

set(fig,'papersize',[6.0 5.5], 'Position',[200 200 900 350]);
print(fig,'plots/equivalent_stiffness_and_beta_different_q','-dpng','-r1000')


%% Equivalent stiffness and damping for diferent non-linearities
str0 = 'duffing';
str1 = 'omegaeq_';
q = 0.75;
str2 = sprintf('fractional_%.2f', q);
str3 = sprintf('mcssamples_%d', 14000);

markers = ["--", "-.", "-"];
letters = ["a", "b", "c"];

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
    
        fig = figure(5);
        for dof = 1:1:ndof
            subplot(2,ndof, dof); 
            hold on
            plot(time,omega_eq_2(dof,:), marker, 'linewidth',2)
            xlabel('Time (s)')
            ylabel('$\omega^2_{e} $ (N/m)')
            xlim([0 4])
            ylim([0 1500])
            xticks([0 1 2 3 4])
            aux = sprintf("%s) DOF %d", letters(dof), dof);
            title(aux)
            grid
            legend('$\epsilon$ = 0.50', '$\epsilon$ = 1.00', '$\epsilon$ = 2.00')
        end
    end
end

str1 = 'betaeq_';
letters = ["d", "e", "f"];

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
    
        figure(5);
        for dof = 1:1:ndof
            subplot(2,ndof, dof+ndof); 
            hold on
            plot(time,beta_eq(dof,:), marker, 'linewidth',2)
            xlim([0 4])
            ylim([0 35])
            xticks([0 1 2 3 4])
            grid
            xlabel('Time (s)')
            ylabel('$\beta_{e} $ (Ns/m)')
            aux = sprintf("%s) DOF %d", letters(dof), dof);
            title(aux)
            legend('$\epsilon$ = 0.50', '$\epsilon$ = 1.00', '$\epsilon$ = 2.00', 'location', 'southeast')
        end
    end
end

set(fig,'papersize',[6.0 5.5], 'Position',[200 200 900 350]);
print(fig,'plots/equivalent_stiffines_and_damping_different_epi','-dpng','-r1000')

%% displacement same q and same ep
load('data/displacement_variance_oscillator_duffing_ndof_3_fractional_0.75_dt_0.0010_mcssamples_14000_damping_40.00_stiffness_400.00_barrier_0.25_powerspectrum_eps_S0_0.20_duffingparameter_epx_1.00.mat')
ndof = 3;

fig = figure(6);
for i=1:ndof
    subplot(1,ndof,i); 
    hold on
    plot(time, sqrt(varx_sl(i,:)),'k-','linewidth',2)
    plot(time, sqrt(c(i,:))','r--','linewidth',2)
    plot(time_out, sqrt(varx_mcs(i,:)),'b:','linewidth',2)
    legend('SL', 'SA', 'MCS', 'location', 'southeast')
    xlabel('Time (s)')
    ylabel('$\sqrt{Var[x(t)]} (m)$')
    xlim([0 4])
    xticks([0 1 2 3 4])
    yticks([0 0.003 0.006 0.009 0.012 0.015 0.018])
    ylim([0 0.018])
    aux = sprintf("DOF %d", i);
    grid(1);
    title(aux)
end

set(fig,'papersize',[6.0 5.5], 'Position',[200 200 900 350]);
print(fig,'plots/displacement_std_same_q_same_epx','-dpng','-r1000')


%% Plot amplitude PDF
lam = 0.25;
str1 = 'pdfs_';
epx = 1;
ndof = 3;
str2 = sprintf('epx_%.2f', epx);
q = 0.75;

letters = ["a" "b" "c" "d" "e" "f"];

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
                fig = figure(3);
                colormap jet
                subplot(2,3,j)
                hold on
                surf(time_out,av,pr(:,:,j));
                plot3(time_out,ones(size(time_out))*barrier(j),ha,'r','linewidth',2)
                clim([0, 300])
                ylim([0, 0.04])
                xlim([0 4])
                xticks([0 1 2 3 4])
                shading interp
                xlabel('Time (s)')
                ylabel('MCS Amplitude (m)')
                aux = sprintf('%s) DOF %d', letters(j), j);
                title(aux)
        
                subplot(2,3,j+ndof)
                hold on
                surf(time_out,av,pa(:,:,j));
                plot3(time_out,ones(size(time_out))*barrier(j),ha,'r','linewidth',2)
                clim([0, 300])
                ylim([0, 0.04])
                xlim([0 4])
                xticks([0 1 2 3 4])
                shading interp
                xlabel('Time (s)')
                ylabel('Approx. Amplitude (m)')
                aux = sprintf('%s) DOF: %d', letters(j+ndof), j);
                title(aux)
            end
            break
        end
    end
end

set(fig,'papersize',[6.0 5.5], 'Position',[200 200 900 350]);
print(fig,'plots/amplitude_pdf','-dpng','-r1000')

%% survival probability
str1 = 'firsttimepassage_';

for i = 1:1:numel(files)
    if (contains(files(i).name, str1) && contains(files(i).name, str2))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(7));

        if (vec(:).q == 0.75 || vec(:).q == 0.50)
            P = vec.P;
            time = vec.time;
            fpp = vec.fpp;
            tfp = vec.tfp;
            files(i).name
    
            fig = figure(4);
            for k = 1:ndof
                if (vec(:).q == 0.75)
                    subplot(2,ndof,k);
                else
                    subplot(2,ndof,k+ndof);
                end
                hold on
                plot(time, P(k,:)','b','linewidth',2);
                plot(tfp, fpp,'r--','linewidth',2);          
                legend('Analytical','MCS')
                aux = sprintf('$q = %.2f$; DOF: %d', vec(:).q, k);
                title(aux)
                xlabel('Time (s)')
                ylabel('Survival propability')
                xlim([0 4])
                ylim([0 1])
                xticks([0 1 2 3 4])
                grid(1)
            end
        end
    end
end

set(fig,'papersize',[6.0 5.5], 'Position',[200 200 900 350]);
print(fig,'plots/survival_prop','-dpng','-r1000')