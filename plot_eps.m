clc
clear
close all

freq_max = 100;
t_max = 30;

freq = 0:0.1:freq_max;
time = 0:0.1:t_max;

eps = zeros(numel(freq), numel(time));
for i=1:1:numel(freq)
    for j = 1:1:numel(time)
        eps(i,j) = evolutionary_power_spectrum(freq(i), time(j));
    end
end

fig = figure(1);
subplot(2, 1, 1)
surf(time, freq, eps)
xlabel('Time (s)', 'FontSize', 12, 'Interpreter','latex')
ylabel('Frequency (rad/s)','FontSize', 12, 'Interpreter','latex')
zlabel('EPS','FontSize', 12, 'Interpreter','latex')
colormap jet
shading interp
view([69 16])

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
set(gca,'TickLabelInterpreter','latex')

subplot(2, 1, 2)
surf(time, freq, eps)
xlabel('Time (s)', 'FontSize', 12, 'Interpreter','latex')
ylabel('Frequency (rad/s)','FontSize', 12, 'Interpreter','latex')
colormap jet
shading interp
view([0 90])

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
set(gca,'TickLabelInterpreter','latex')

%saveas(fig, 'plots/eps', 'pdf')
set(fig,'papersize',[6.75 5.0]);
print(fig,'plots/eps','-dpdf')