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

set(groot,'defaultAxesFontSize', 16)

%% Equivalent stiffness and damping
a = 0.50;
barrier = 0.25;
xy = 0.001;
str0 = sprintf('xy_%.2f', xy);
str1 = 'omegaeq_';
str2 = sprintf('_a_%.2f', a);
str3 = sprintf('barrier_%.2f', barrier);

markers = ["-", "-."];
letters = ["a", "b", "b"];

for i = 1:1:numel(files)
    if(contains(files(i).name, str0) && contains(files(i).name, str1) && contains(files(i).name, str2) && contains(files(i).name, str3))
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
                xlabel('Time')
                ylabel('$\omega^2_{e}$')
                aux = sprintf("%s) DOF %d", letters(dof), dof);
                title(aux, 'fontsize', 18)
                grid(1);
                xlim([0 4])
                ylim([0 1600])
                xticks([0 1 2 3 4])
                legend('q = 0.5', 'q = 0.75')
            end
        end
    end
end

str1 = 'betaeq_';
letters = ["d", "e", "f"];

for i = 1:1:numel(files)
    if(contains(files(i).name, str0) && contains(files(i).name, str1) && contains(files(i).name, str2) && contains(files(i).name, str3))
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
                xlabel('Time')
                ylabel('$\beta_{e}$')
                aux = sprintf("%s) DOF %d", letters(dof), dof);
                title(aux, 'fontsize', 18)
                grid(1);
                xlim([0 4])
                ylim([0 45])
                xticks([0 1 2 3 4])
                legend('q = 0.5', 'q = 0.75', 'location', 'north')
            end
        end
    end
end

set(fig,'papersize',[6.0 5.5], 'Position',[200 200 900 350]);
print(fig,'plots/equivalent_stiffness_and_beta_different_q_bw','-dpng','-r1000')

%% Equivalent stiffness and damping for diferent non-linearities
q = 0.75;
str0 = sprintf('xy_%.2f', xy);
str1 = 'omegaeq_';
str2 = sprintf('fractional_%.2f', q);
str3 = sprintf('mcssamples_%d', 14000);
str4 = sprintf('barrier_%.2f', barrier);

markers = ["--", "-.", "-"];
letters = ["a", "b", "c"];

for i = 1:1:numel(files)
    if(contains(files(i).name, str0) && contains(files(i).name, str1) ...
            && contains(files(i).name, str2) && contains(files(i).name, str3) && contains(files(i).name, str4))
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
            xlabel('Time')
            ylabel('$\omega^2_{e}$')
            xlim([0 4])
            ylim([0 1500])
            xticks([0 1 2 3 4])
            aux = sprintf("%s) DOF %d", letters(dof), dof);
            title(aux, 'fontsize', 18)
            grid
            legend('$\alpha$ = 0.30', '$\alpha$ = 0.50', '$\alpha$ = 0.70')
        end
    end
end

str1 = 'betaeq_';
letters = ["d", "e", "f"];

for i = 1:1:numel(files)
    if(contains(files(i).name, str0) && contains(files(i).name, str1) ...
            && contains(files(i).name, str2) && contains(files(i).name, str3) && contains(files(i).name, str4))
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
            ylim([0 60])
            xticks([0 1 2 3 4])
            grid
            xlabel('Time (s)')
            ylabel('$\beta_{e} $ (Ns/m)')
            aux = sprintf("%s) DOF %d", letters(dof), dof);
            title(aux, 'fontsize', 18)
            legend('$\alpha$ = 0.30', '$\alpha$ = 0.50', '$\alpha$ = 0.70', 'location', 'north')
        end
    end
end

set(fig,'papersize',[6.0 5.5], 'Position',[200 200 900 350]);
print(fig,'plots/equivalent_stiffines_and_damping_different_epi_bw','-dpng','-r1000')

%% Plot amplitude PDF
lam = 0.25;
a = 0.50;
ndof = 3;
q = 0.75;
str1 = 'pdfs_';
str2 = sprintf('_a_%.2f', a);
str3 = sprintf('barrier_%.2f', lam);
str4 = sprintf('xy_%.2f', xy);

letters = ["a" "b" "c" "d" "e" "f"];

for i = 1:1:numel(files)
    if (contains(files(i).name, 'displacement_variance') && contains(files(i).name, str2) && contains(files(i).name, str3) && contains(files(i).name, str4))
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
    
    if(contains(files(i).name, str1) && contains(files(i).name, str2) && contains(files(i).name, str3) && contains(files(i).name, str4))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(7));

        if (str2double(str_split(7)) == q)
            av = vec.av;
            pa = vec.pa;
            pr = vec.pr;
            time_out = vec.time_out;

            ha = ones(size(time_out))*1000;
            fig = figure(3);

            for j = 1:1:ndof
                colormap jet
                subplot(2,3,j)
                hold on
                surf(time_out,av,pr(:,:,j));
                plot3(time_out,ones(size(time_out))*barrier(j),ha,'r','linewidth',2)
                clim([0, 300])
                ylim([0, 0.03])
                xlim([0 4])
                xticks([0 1 2 3 4])
                shading interp
                xlabel('Time')
                ylabel('Amplitude')
                aux = sprintf('%s) DOF %d', letters(j), j);
                title(aux, 'fontsize', 18)
        
                subplot(2,3,j+ndof)
                hold on
                surf(time_out,av,pa(:,:,j));
                plot3(time_out,ones(size(time_out))*barrier(j),ha,'r','linewidth',2)
                clim([0, 300])
                ylim([0, 0.03])
                xlim([0 4])
                xticks([0 1 2 3 4])
                shading interp
                xlabel('Time')
                ylabel('Amplitude')
                aux = sprintf('%s) DOF: %d', letters(j+ndof), j);
                title(aux, 'fontsize', 18)
            end
            break
        end    
    end
end

set(fig,'papersize',[6.0 5.5], 'Position',[200 200 900 350]);
print(fig,'plots/amplitude_pdf_bw','-dpng','-r1000')

%% survival probability
barrier = 0.25;
str1 = 'firsttimepassage_';
str2 = sprintf('_a_%.2f', a);
str3 = sprintf('barrier_%.2f', barrier);
str4 = sprintf('xy_%.2f', xy);

for i = 1:1:numel(files)
    if (contains(files(i).name, str1) && contains(files(i).name, str2) && contains(files(i).name, str3) && contains(files(i).name, str4))
        vec = load(strcat('data/', files(i).name));
        str_split = strsplit(files(i).name,"_");
        vec(:).q = str2double(str_split(7));

        if (vec(:).q == 0.75 || vec(:).q == 0.50)
            P = vec.P;
            time = vec.time;
            fpp = vec.fpp;
            tfp = vec.tfp;
    
            fig = figure(4);
            for k = 1:ndof
                if (vec(:).q == 0.75)
                    subplot(2,ndof,k);
                    aux = sprintf('%s) $q = %.2f$; DOF: %d', letters(k), vec(:).q, k);
                else
                    subplot(2,ndof,k+ndof);
                    aux = sprintf('%s) $q = %.2f$; DOF: %d', letters(k+ndof), vec(:).q, k);
                end
                hold on
                plot(time, P(k,:)','b','linewidth',2);
                plot(tfp, fpp,'r--','linewidth',2);          
                legend('Analytical','MCS')
                title(aux, 'fontsize', 18)
                xlabel('Time')
                ylabel('Survival propability')
                xlim([0 2])
                ylim([0 1])
                xticks([0 0.5 1 1.5 2])
                grid(1)
            end
        end
    end
end

set(fig,'papersize',[6.0 5.5], 'Position',[200 200 900 350]);
print(fig,'plots/survival_prop_bw','-dpng','-r1000')