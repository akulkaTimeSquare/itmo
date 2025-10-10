[y, Fs] = audioread('ac.mp3');
y = y(:,1);
y = y(:).';
factor = 2;         % или 20, 50 — по ситуации
y = y(1:factor:end); % уменьшаем длину сигнала
Fs = Fs / factor;
t = (0:length(y)-1)/Fs;

f = figure();
%plot(t, y, 'b', 'LineWidth', 5);
%xlabel('Время', 'FontSize', 20); ylabel('Значения', "FontSize", 20);
%title('Исходный сигнал', 'FontSize', 25);
%grid on;

V = 1000;              % Максимальная интересующая частота (Гц)
dv = 2;                % Шаг по частоте
v = 0:dv:V;           % Частоты от -V до V

Y = zeros(size(v));
for k = 1:length(v)
    Y(k) = trapz(t, y .* exp(-1i*2*pi*v(k)*t));
end

[peaks, locs] = findpeaks(abs(Y), v, 'MinPeakHeight', 0.125*max(abs(Y)));

peaks
locs

magY = abs(Y);
plot(v, magY, 'r', 'LineWidth', 3); hold on;
plot(locs, peaks, 'ko');
xlabel('Частота', "FontSize", 20); ylabel('Амплитуда', 'FontSize', 20);
title('Модуль Фурье-образа', 'FontSize', 25);
grid on;

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [100, 100, 1800, 1000]);
set(gca, 'FontSize', 22, 'GridLineWidth', 3, 'LineWidth', 2);
exportgraphics(f, "images\16.jpg");

