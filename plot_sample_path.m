clc
clear
close all

T = 30;
freq = 10;
n_paths = 2;
time = 0:0.01:T;
m=800;
frequency_max = 50;
dw = frequency_max/m;

W = zeros(numel(time),n_paths);

for path = 1:1:n_paths
    for ii=0:m-1
        eps_fun = evolutionary_power_spectrum(ii*dw,time);
        W(:,path) = W(:,path) + (((4*dw.*eps_fun).^0.5).*cos((ii*dw).*time - (2*pi)*rand))'; 
    end
end

subplot(2,1,1);
hAx=gca;
hAx.FontName='Times New Roman';
hAx.FontSize=15;
plot(time, W(:,1), 'LineWidth', 1.5)
title('Realizations of a stochastic process', 'FontSize', 16)
xlabel('Time', 'FontSize', 14)
ylabel('w(t)', 'FontSize', 14)

subplot(2,1,2);
plot(time, W(:,2), 'LineWidth', 1.5)
xlabel('Time', 'FontSize', 14)
ylabel('w(t)', 'FontSize', 14)