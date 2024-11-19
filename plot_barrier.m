clc
clear
close all

n = 100;
vec = zeros(1, n);
barrier = 0.9;
time = 0:1:n-1;

a = -1;
b = 1;

for dw = 1:1:n
    vec(1, dw) = a + (b-a).*rand;
end

passage = vec(vec(1,:) >= barrier | vec(1,:) <= -barrier);
first_time = time(vec(1,:) >= barrier | vec(1,:) <= -barrier);

plot(time, vec,'LineWidth', 1.5)
hAx=gca;
hAx.FontName='Times New Roman';
hAx.FontSize=15;
hold on
plot(time, barrier*ones(1, n), 'r', 'LineWidth', 1)
plot(time, -barrier*ones(1, n), 'r', 'LineWidth', 1)
scatter(first_time, passage,100,'r', marker='.')
title('Realization of a stochastic process', 'FontSize', 20)
xlabel('Time', 'FontSize', 18)
ylabel('Amplitude', 'FontSize', 18)
legend('Process', 'Barrier')
ylim([-1.1 1.1])
