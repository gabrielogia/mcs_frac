clc
clear
close all

freq_max = 50;
t_max = 20;

freq = 0:0.01:freq_max;
time = 0:0.01:t_max;

eps = zeros(numel(freq), numel(time));
for i=1:1:numel(freq)
    for j = 1:1:numel(time)
        eps(i,j) = evolutionary_power_spectrum(freq(i), time(j));
    end
end

surf(freq, time, eps')
xlabel('Time (s)', 'FontSize', 14)
ylabel('Frequency (rad/s)','FontSize', 14)
zlabel('EPS','FontSize', 14)
title('Nonseparable excitation evolutionary power spectrum', 'FontSize', 16)
colormap jet
shading interp