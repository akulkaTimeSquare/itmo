% Параметры
T = 4;            % Общая длина временного окна
dt = 0.0075;      % Шаг дискретизации
t = -T/2:dt:T/2;

nt = -T/2:0.001:T/2;
np = zeros(size(nt));
np(abs(nt) <= 1/2) = 1;

% Исходная функция
Pi_t = double(abs(t) <= 0.5);
N = length(Pi_t);
% Прямое БПФ
Pi_f = fftshift(fft(Pi_t)) / sqrt(N);

% Частотная ось
dnu = 1/T;
V = 1/dt;
nu = -V/2:dnu:V/2;

nnu = -V/2:0.01:V/2;
hatp = sinc(nnu);

% Обратное БПФ
Pi_t_rec = ifft(ifftshift(Pi_f)) * sqrt(N);
f = figure();

plot(t, Pi_t_rec, 'b', 'LineWidth', 7); hold on;
plot(nt, np, 'k--', LineWidth=6);
xlabel('Время t', 'FontSize', 20); ylabel('Значения', "FontSize", 20);
title('Восстановленная прямоугольная функция П(t)', 'FontSize', 25);
legend("Восстановленная функция", "Исходная функция", FontSize = 17);
ylim([-0.1 1.1]);
xlim([-1 1]);
grid on;

%plot(nu, Pi_f, 'r', 'LineWidth', 9); hold on;
%plot(nnu, hatp, 'k', 'LineWidth', 4);
%%xlabel('Частота', "FontSize", 20); ylabel('Амплитуда', 'FontSize', 20);
%title('Фурье-образ П(t) (быстрое преобразование Фурье)', 'FontSize', 25);
%legend("Вычисленный Фурье-образ", "Аналитический Фурье-образ", FontSize = 17)
%xlim([-5 5]);
%grid on;

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [100, 100, 1800, 1000]);
set(gca, 'FontSize', 22, 'GridLineWidth', 3, 'LineWidth', 2);
%exportgraphics(f, "images\35.jpg")