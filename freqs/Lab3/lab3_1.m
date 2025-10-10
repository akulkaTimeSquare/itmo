[y, f] = audioread('MUHA.wav'); % y — сигнал, f — частота дискретизации
%sound(y, f);                    % прослушать аудио

dt = 1/f;
T = length(y)*dt;
t = 0:dt:T-dt;

Y = fftshift(fft(y));
N = length(Y);
%nu = (0:N-1)*(f/N); % ось частот
magY = abs(Y);          % модуль спектра
nu = (-N/2:N/2-1)*(f/N);

hatnu = 5;
n1 = 300;
n2 = 5000;
Y_filtered = Y;
Y_filtered(nu > -n1 & nu < n1) = 0;
Y_filtered(nu < -n2 | nu > n2) = 0;

y_filtered = real(ifft(ifftshift(Y_filtered)));
sound(y_filtered, f);

f = figure();
%plot(nu, magY, 'r', 'LineWidth', 6); hold on;
%plot(nu, abs(Y_filtered), 'g', 'LineWidth', 4);
%legend('Модуль исходного', 'Модуль отфильтрованного', 'FontSize', 21);
%xlabel('Частота', "FontSize", 20); ylabel('Амплитуда', 'FontSize', 20);
%title('Модули Фурье-образов', 'FontSize', 25);
%xlim([-6500, 6500])
%grid on;

plot(t, y, 'b', 'LineWidth', 5); hold on;
plot(t, y_filtered, 'g', 'LineWidth', 8);
xlabel('Время', 'FontSize', 20); ylabel('Значения', "FontSize", 20);
title('Исходный сигнал', 'FontSize', 25);
legend('Зашумленный сигнал', 'Отфильтрованный', "FontSize", 22);
grid on;

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [100, 100, 1800, 1000]);
set(gca, 'FontSize', 22, 'GridLineWidth', 3, 'LineWidth', 2);
%exportgraphics(f, "images\40.jpg")