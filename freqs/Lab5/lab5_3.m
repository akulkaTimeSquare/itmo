T = 1.15;
dt = 0.05;
N = T/dt;
t = linspace(-T/2, T/2, N);
%p = exp(t.*double(abs(t) <= 0.5)); % прямоугольная функция
p = double(abs(t) <= 0.5);

nt = -T/2:0.0001:T/2;
np = zeros(size(nt));
np(abs(nt) <= 1/2) = 1;

% FFT и спектр
m = 0:N-1;
c = dt * (-1).^m;
V = 1/dt;
dnu = 1/T;
nu = linspace(-V/2, V/2-dnu, N);

F_fft = fftshift(c .* fft(p));

% Аналитический спектр
nnu = -V/2:0.001:V/2;
hatp = sinc(nnu);

% Обратное преобразование Фурье
f_recovered = real(ifft(ifftshift(F_fft) ./ c));


F_plain = fftshift(fft(p)) / sqrt(N); % Без учета c_m
f_recovered_with_ifft = ifft(ifftshift(F_plain)) * sqrt(N);


hdnu = 1/T;
hV = 1/dt;
hnu = nu;

Pi_nu = zeros(size(hnu));
for k = 1:length(hnu)
    Pi_nu(k) = trapz(t, p .* exp(-1j*2*pi*hnu(k)*t));
end

% Обратное преобразование Фурье
Pi_t_rec = zeros(size(t));
for k = 1:length(t)
    Pi_t_rec(k) = trapz(hnu, Pi_nu .* exp(1j*2*pi*hnu*t(k)));
end

% Нормировка (по аналогии с непрерывным интегралом)
Pi_t_rec = real(Pi_t_rec);


f = figure();
%plot(t, f_recovered, 'r', 'LineWidth', 16); hold on;
%plot(t, f_recovered_with_ifft, 'g', 'LineWidth', 7); hold on;
%plot(t, Pi_t_rec, 'b', 'LineWidth', 7); hold on;
%plot(nt, np, 'k--', LineWidth=6);
%xlabel('Время t', 'FontSize', 20); ylabel('Значения', "FontSize", 20);
%title('Сравнение методов при восстановлении прямоугольной функции', 'FontSize', 25);
%legend("Улучшенный fft восстановление", "fft восстановление", ...
%    "Численное интегрирование восстановление", "Исходная функция", FontSize = 17);
%ylim([-0.1 1.1]);
%xlim([-T/2 T/2]);
%grid on;


plot(nu, F_fft, 'r', 'LineWidth', 5); hold on;
plot(nu, F_plain, 'g', 'LineWidth', 3);
plot(hnu, Pi_nu, 'b', 'LineWidth', 6);
plot(nnu, hatp, 'k', 'LineWidth', 4);
xlabel('Частота', "FontSize", 20); ylabel('Амплитуда', 'FontSize', 20);
title('Фурье-образ П(t) (улучшенное быстрое преобразование Фурье)', 'FontSize', 25);
legend("Улучшенный fft Фурье-образ", "fft Фурье-образ", ...
    "Численное интегрирование Фурье-образ", "Анатилический Фурье-образ", FontSize = 17);
xlim([-10/2+dnu 10/2-dnu]);
ylim([-0.25 1.1])
grid on;

sum(abs(abs(sinc(nu)) - abs(F_fft)))
sum(abs(abs(sinc(hnu)) - abs(Pi_nu)))

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [100, 100, 1800, 1000]);
set(gca, 'FontSize', 22, 'GridLineWidth', 3, 'LineWidth', 2);
%exportgraphics(f, "images\56.jpg")