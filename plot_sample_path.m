clc
clear
close all

T = 1;
freq = 10;
n_paths = 2;
time = 0:0.01:T;
vec = zeros(n_paths, numel(time));

a = -1;
b = 1;

for path = 1:1:n_paths
    for dt = 1:1:numel(time)
        vec(path, dt) = (a + (b-a).*rand)*evolutionary_power_spectrum(freq, dt);
    end
end

subplot(2,1,1);
hAx=gca;
hAx.FontName='Times New Roman';
hAx.FontSize=15;
plot(time, vec(1,:), 'LineWidth', 1.5)
title('Realizations of a stochastic process', 'FontSize', 16)
xlabel('Time', 'FontSize', 14)
ylabel('Amplitude', 'FontSize', 14)

subplot(2,1,2);
plot(time, vec(2,:), 'LineWidth', 1.5)
xlabel('Time', 'FontSize', 14)
ylabel('Amplitude', 'FontSize', 14)