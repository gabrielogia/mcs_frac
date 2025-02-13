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
            vec(:).q
            for dof = 1:1:ndof
                subplot(2,ndof, dof); 
                hold on
                plot(time,omega_eq_2(dof,:), marker, 'linewidth',2)
                xlabel('Time')
                ylabel('$\omega^2_{e} (t)$')
                aux = sprintf("%s) DOF %d", letters(dof), dof);
                title(aux,'FontSize',18)
                grid("on");
                xlim([0 4])
                ylim([0 1500])
                xticks([0 1 2 3 4])
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
                xlabel('Time')
                ylabel('$\beta_{e} (t)$')
                aux = sprintf("%s) DOF %d", letters(dof), dof);
                title(aux, 'FontSize',18)
                grid("on");
                xlim([0 4])
                xticks([0 1 2 3 4])
                ylim([0 50])
                yticks([0 25 50])
                legend('q = 0.5', 'q = 0.75', 'location', 'north')
            end
        end
    end
end

set(fig,'papersize',[6.0 5.5], 'Position',[200 200 900 350]);
print(fig,'plots/equivalent_stiffness_and_beta_different_q','-dpng','-r1000')

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
                ylim([0, 0.03])
                xlim([0 4])
                xticks([0 1 2 3 4])
                shading interp
                xlabel('Time')
                ylabel('Amplitude')
                aux = sprintf('%s) DOF %d', letters(j), j);
                title(aux)
        
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
                aux = sprintf('q$ = %.2f$; DOF: %d', vec(:).q, k);
                title(aux)
                xlabel('Time')
                ylabel('Survival propability')
                xlim([0 2])
                ylim([0 1])
                xticks([0 1 2 3 4])
                grid(1)
            end
        end
    end
end

set(fig,'papersize',[6.0 5.5], 'Position',[200 200 900 350]);
print(fig,'plots/survival_prop','-dpng','-r1000')