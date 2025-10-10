% === Параметры для анализа ===
T = 4;
dt = 0.05;
N = T/dt;
t = linspace(-T/2, T/2-dt, N);
func = double(abs(t) <= 0.5); % прямоугольная функция

nt = -T/2:0.0001:T/2;
np = zeros(size(nt));
np(abs(nt) <= 1/2) = 1;

% FFT и спектр
m = linspace(-N/2, N/2-1, N);
c = dt * (-1).^m;
V = 1/dt;
dnu = 1/T;
nu = m / T;

F_fft = fftshift(c .* fft(func));

% Аналитический спектр
nnu = -V/2:0.001:V/2;
hatp = sinc(nnu);

% Обратное преобразование Фурье
f_recovered = real(ifft(ifftshift(F_fft) ./ c));


f = figure();

plot(t, f_recovered, 'b', 'LineWidth', 7); hold on;
plot(nt, np, 'k--', LineWidth=6);
xlabel('Время t', 'FontSize', 20); ylabel('Значения', "FontSize", 20);
title('Восстановленная прямоугольная функция П(t)', 'FontSize', 25);
legend("Восстановленная функция", "Исходная функция", FontSize = 17);
ylim([-0.1 1.1]);
xlim([-2/2 2/2]);
grid on;

%plot(nu, F_fft, 'r', 'LineWidth', 9); hold on;
%plot(nnu, hatp, 'k--', 'LineWidth', 6);
%xlabel('Частота', "FontSize", 20); ylabel('Амплитуда', 'FontSize', 20);
%title('Фурье-образ П(t) (улучшенное быстрое преобразование Фурье)', 'FontSize', 25);
%legend("Вычисленный Фурье-образ", "Аналитический Фурье-образ", FontSize = 17)
%xlim([-V/2+dnu V/2-dnu]);
%grid on;


set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [100, 100, 1800, 1000]);
set(gca, 'FontSize', 22, 'GridLineWidth', 3, 'LineWidth', 2);
exportgraphics(f, "images\35.jpg")