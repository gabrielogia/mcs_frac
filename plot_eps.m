clc
clear
close all

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

set(groot,'defaultAxesFontSize',12)

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
set(fig,'papersize',[7 5.0]);
print(fig,'plots/eps','-dpdf')