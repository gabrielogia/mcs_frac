clc
clear
close all

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

set(groot,'defaultAxesFontSize', 16)

freq_max = 100;
t_max = 30;

freq = 0:0.1:freq_max;
time = 0:0.1:t_max;

S0 = 0.2;

eps = zeros(numel(freq), numel(time));
for i=1:1:numel(freq)
    for j = 1:1:numel(time)
        eps(i,j) = evolutionary_power_spectrum(freq(i), time(j), S0);
    end
end

fig = figure(1);
subplot(1, 2, 1)
surf(time, freq, eps)
xlabel('Time (s)')
ylabel('Frequency (rad/s)')
zlabel('EPS')
colormap jet
shading interp
view([65 11])

subplot(1, 2, 2)
surf(time, freq, eps)
xlabel('Time (s)')
ylabel('Frequency (rad/s)')
colormap jet
shading interp
view([0 90])

set(fig,'papersize',[6.0 5.5], 'Position', [200 200 900 350]);
print(fig,'plots/eps','-dpng','-r1000')