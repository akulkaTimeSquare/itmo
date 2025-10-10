% Параметры для y1(t)
a1 = 1;
a2 = 0.5;
w1 = 2*pi*5; % 5 Гц
w2 = 2*pi*10; % 10 Гц
phi1 = 0;
phi2 = pi/4;
B1 = 20;

% Временной массив для "непрерывной" функции
T = 1e3;
T1 = 5;
dt = 1e-3; % очень маленький шаг для имитации непрерывности
t = -T/2:dt:T/2; % берем чуть больше интервала для хорошей картинки
N = length(t);
m = 0:N-1;
c = dt * (-1).^m;

% Определение функций
y1 = a1*sin(w1*t + phi1) + a2*sin(w2*t + phi2);

y1_fft = fftshift(c .* fft(y1));

f = figure();
%plot(t, y1, 'b', 'LineWidth', 6); hold on;
%xlabel('Время', 'FontSize', 20); ylabel('Значения', "FontSize", 20);
%title('График функции y_1(t)', 'FontSize', 25);
%xlim([-T1/2 T1/2]);
%ylim([-1.75 1.25]);
%grid on;

Fs = 30; % Частота дискретизации в Гц
dt_sample = 1/Fs;
t_sample = -T1/2:dt_sample:T1/2;

% Сэмплированные значения
y1_sample = a1*sin(w1*t_sample + phi1) + a2*sin(w2*t_sample + phi2);

y1_interp = zeros(size(t));
for k = 1:length(t_sample)
    y1_interp = y1_interp + y1_sample(k)*sinc(Fs*(t - t_sample(k)));
end

y1_fft_interp = fftshift(c .* fft(y1_interp));

c1_sample = dt_sample * (-1).^(0:length(t_sample)-1);
y1_sample_fft = fftshift(c1_sample .* fft(y1_sample));

%plot(t, y1_interp, 'g', LineWidth=8); hold on;
%plot(t, y1, 'k', LineWidth=4);
%scatter(t_sample, y1_sample, "filled", MarkerFaceColor="red", MarkerEdgeColor="red", LineWidth=5)
%title('Восстановление функции через теорему Котельникова');
%xlabel('Время');
%ylabel('Значения');
%legend('Интерполяция', 'Оригинальная', "Сэмплированная");
%xlim([-T1/2-0.25 T1/2+0.25]);
%grid on;


V = 1/dt_sample;
dnu = 1/T1;
nu = -1/dt/2:1/T:1/dt/2;
N_sample = T1/dt_sample+1;
nu_sample = linspace(-V/2, V/2, N_sample);

plot(nu, y1_fft_interp, 'g', 'LineWidth', 12); hold on;
plot(nu, y1_fft, 'k', 'LineWidth', 6);
plot(nu_sample, y1_sample_fft, 'r', 'LineWidth', 8);
plot(nu_sample-V, y1_sample_fft, 'r', 'LineWidth', 8, 'HandleVisibility', 'off');
plot(nu_sample+V, y1_sample_fft, 'r', 'LineWidth', 8, 'HandleVisibility', 'off');
xlabel('Частота', "FontSize", 20);
ylabel('Амплитуда', 'FontSize', 20);
title('Сравнение образов Фурье при применении теоремы Котельникова', 'FontSize', 25);
xlim([-B1/2-5 B1/2+5]);
xline([-B1/2, B1/2], '--b', 'B = 20', 'FontSize', 20, 'LineWidth', 4, 'HandleVisibility', 'off');
xline([-V/2, V/2], '--b', 'V = ' + string(V), 'FontSize', 20, 'LineWidth', 4, 'HandleVisibility', 'off');
legend('Фурье-образ интерполяции', 'Фурье-образ оригинала', ...
    "Фурье-образ сэмплированного сигнала", Location='southeast');
grid on;
ylim([-0.5 0.5]);

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [100, 100, 1800, 1000]);
set(gca, 'FontSize', 22, 'GridLineWidth', 3, 'LineWidth', 2);
%exportgraphics(f, "images\67.jpg")