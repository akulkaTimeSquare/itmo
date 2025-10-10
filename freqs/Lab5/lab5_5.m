% Параметры для y2(t)
b = 4*pi; % b определяет ширину основного лепестка у sinc(bt)
B2 = 4;

% Временной массив для "непрерывной" функции
T = 1e2;
T2 = 15; % промежуток для отображения, можно варьировать
dt = 1e-3; % очень маленький шаг для имитации непрерывности
N = T/dt;
t = linspace(-T/2, T/2, N); % берем чуть больше интервала для хорошей картинки
N = length(t);
m = 0:N-1;
c = dt * (-1).^m;

% Определение функций
y2 = sinc(b*t/pi); % В MATLAB sinc(x) = sin(pi*x)/(pi*x), поэтому делим на pi

y2_fft = fftshift(c .* fft(y2));

f = figure();
%plot(t, y2, 'b', 'LineWidth', 7);
%xlabel('Время', 'FontSize', 20); ylabel('Значения', "FontSize", 20);
%title('График функции y_2(t)', 'FontSize', 25);
%xlim([-T2/2 T2/2]);
%ylim([-0.5 1.15]);
%grid on;

Fs = 10; % Частота дискретизации в Гц
dt_sample = 1/Fs;
N_sample = T2/dt_sample;
t_sample = linspace(-T2/2, T2/2, N_sample);

% Сэмплированные значения
y2_sample = sinc(b*t_sample/pi);

y2_interp = zeros(size(t));
for k = 1:length(t_sample)
    y2_interp = y2_interp + y2_sample(k)*sinc(Fs*(t - t_sample(k)));
end

y2_fft_interp = fftshift(c .* fft(y2_interp));

c2_sample = dt_sample * (-1).^(0:length(t_sample)-1);
y2_sample_fft = fftshift(c2_sample .* fft(y2_sample));

%plot(t, y2_interp, 'g', LineWidth=8); hold on;
%plot(t, y2, 'k', LineWidth=4);
%scatter(t_sample, y2_sample, "filled", MarkerFaceColor="red", MarkerEdgeColor="red", LineWidth=1)
%title('Восстановление функции через теорему Котельникова');
%xlabel('Время');
%ylabel('Значения');
%legend('Интерполяция', 'Оригинальная', "Сэмплированная");
%xlim([-T2/2-0.25 T2/2+0.25]);
%ylim([-0.5 1.15]);
%grid on;


V = 1/dt_sample;
dnu = 1/T2;
nu = linspace(-1/dt/2, 1/dt/2, N);
nu_sample = linspace(-V/2, V/2, N_sample);

plot(nu, y2_fft_interp, 'g', 'LineWidth', 12); hold on;
plot(nu, y2_fft, 'k', 'LineWidth', 4);
plot(nu_sample, y2_sample_fft, 'r', 'LineWidth', 8);
plot(nu_sample-V, y2_sample_fft, 'r', 'LineWidth', 8, 'HandleVisibility', 'off');
plot(nu_sample+V, y2_sample_fft, 'r', 'LineWidth', 8, 'HandleVisibility', 'off');
xlabel('Частота', "FontSize", 20);
ylabel('Амплитуда', 'FontSize', 20);
title('Сравнение образов Фурье при применении теоремы Котельникова', 'FontSize', 25);
xlim([-V/2-1 V/2+1]);
ylim([-0.5 0.5]);
xline([-B2/2, B2/2], '--b', 'B = 4', 'FontSize', 20, 'LineWidth', 4, 'HandleVisibility', 'off');
xline([-V/2, V/2], '--b', 'V = ' + string(V), 'FontSize', 20, 'LineWidth', 4, 'HandleVisibility', 'off');
legend('Фурье-образ интерполяции', 'Фурье-образ оригинала', ...
    "Фурье-образ сэмплированного сигнала", Location='southeast');
grid on;

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [100, 100, 1800, 1000]);
set(gca, 'FontSize', 22, 'GridLineWidth', 3, 'LineWidth', 2);
exportgraphics(f, "images\74.jpg")