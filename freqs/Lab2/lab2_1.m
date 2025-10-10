T = 100;
N = 16384;
t = linspace(-T, T, N);
dt = t(2) - t(1);
freq = (-N/2:N/2-1)/(N*dt);
omega = 2*pi*freq;         


% Обновлённые значения параметра c ≠ 0
c_values = [-3, -1, 2, 3];
colors = ['b', 'r', 'g', 'cyan']; % Цвета для графиков
widths = [13, 9, 4, 3];

%%
f = figure()
hold on; grid on;
% Построение спектров
for i = 1:length(c_values)
    c = c_values(i);
    g = exp(-(t + c).^2);     
    plot(t, g, colors(i), 'LineWidth', widths(i), 'DisplayName', ...
         ['g(t) при c = ' num2str(c)]);
end
xlabel('t', 'FontSize', 20); ylabel('g(t)', 'FontSize', 20);
xlim([-5 6]);
ylim([-0.05 1.15]);
title('Смещенные функции g(t)', 'FontSize', 25);
legend show;
hold off;
%%
f = figure()
hold on; grid on;
% Построение спектров
for i = 1:length(c_values)
    c = c_values(i);
    hatg = exp(1i*omega*c).*exp(-omega.^2/4)/sqrt(2);
    rhatg = real(hatg);     
    plot(omega, rhatg, colors(i), 'LineWidth', widths(i), 'DisplayName', ...
         ['c = ' num2str(c)]);
end
xlabel('t', 'FontSize', 20); ylabel('g(t)', 'FontSize', 20);
xlim([-6 6]);
ylim([-0.7 0.8]);
title('Вещественная часть образа g(t)', 'FontSize', 25);
legend show;
hold off;
%%
f = figure()
hold on; grid on;
% Построение спектров
for i = 1:length(c_values)
    c = c_values(i);
    hatg = exp(1i*omega*c).*exp(-omega.^2/4)/sqrt(2);
    imhatg = imag(hatg);     
    plot(omega, imhatg, colors(i), 'LineWidth', widths(i), 'DisplayName', ...
         ['c = ' num2str(c)]);
end
xlabel('t', 'FontSize', 20); ylabel('g(t)', 'FontSize', 20);
xlim([-6 6]);
ylim([-0.7 0.8]);
title('Мнимая часть образа g(t)', 'FontSize', 25);
legend show;
hold off;
%%
f = figure()
hold on; grid on;
% Построение спектров
for i = 1:length(c_values)
    c = c_values(i);
    hatg = exp(1i*omega*c).*exp(-omega.^2/4)/sqrt(2);
    abshatg = abs(hatg);     
    plot(omega, abshatg, colors(i), 'LineWidth', widths(i), 'DisplayName', ...
         ['c = ' num2str(c)]);
end
xlabel('t', 'FontSize', 20); ylabel('g(t)', 'FontSize', 20);
xlim([-6 6]);
ylim([-0.15 0.8]);
title('Модуль образа g(t)', 'FontSize', 25);
legend show;
hold off;
%%
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [100, 100, 1800, 1000]);
set(gca, 'FontSize', 22, 'GridLineWidth', 3, 'LineWidth', 2);
exportgraphics(f, "images\14.jpg")