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

q_vec = [0.5 0.75 1];

str1 = 'powerspectrum_wn';
str2 = 'displacement_';
plot_formulation_diff(files, str1, str2, q_vec, 1)

str1 = 'powerspectrum_eps';
plot_formulation_diff(files, str1, str2, q_vec, 2)

function plot_formulation_diff(files, str1, str2, q_vec, n_fig)
    for i = 1:1:numel(files)
        if(contains(files(i).name, str1) && contains(files(i).name, str2))
            vec = load(strcat('data/', files(i).name));
            str_split = strsplit(files(i).name,"_");
    
            q = str2double(str_split(8));
            formulation = char(str_split(20));
    
            time = vec(:).time;
            time_out = vec(:).time_out;
            c = vec(:).c;
            varx_mcs = vec(:).varx_mcs;
            varx_sl = vec(:).varx_sl;
    
            for j = 1:numel(q_vec)
                fig = figure(n_fig);
                if q == q_vec(j)
                    if formulation == "optimization" 
                        subplot(2,3,j)
                        plot(time, sqrt(c), time, sqrt(varx_sl), time_out, sqrt(varx_mcs), 'LineWidth', 2)
                        aux = sprintf('q = %.2f', q);
                        title(aux)
                        ylabel('$\sigma(y)$ (m)')
                        xlabel('Time (s)')
                        legend('Integral', 'SL', 'MCS', 'Location', 'southeast')
                    else
                        subplot(2,3,j+3)
                        plot(time, sqrt(c), time, sqrt(varx_sl), time_out, sqrt(varx_mcs), 'LineWidth', 2)
                        aux = sprintf('q = %.2f', q);
                        ylabel('$\sigma(y)$ (m)')
                        xlabel('Time (s)')
                        title(aux)
                        legend('WN Approximation', 'SL', 'MCS', 'Location', 'southeast')
                    end
                end
            end
        end
    end
end