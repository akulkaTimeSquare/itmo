a = 1;
t1 = -3;
t2 = 3;
T = 100;
dt = 0.01;
t = -T/2:dt:T/2;
N = length(t);
g = zeros(size(t));
g(t >= t1 & t <= t2) = a;

c = 0.1;
b = 0.35;
d = 3;
u = g + b*(rand(size(t))-0.5) + c*sin(d*t);
%%
f = figure();
plot(t, u, 'b', 'LineWidth', 5); hold on;
plot(t, g, 'k--', 'LineWidth', 8);
xlabel('Время', 'FontSize', 20); ylabel('Значения', "FontSize", 20);
title('Исходный и зашумленный сигналы', 'FontSize', 25);
legend('Зашумленный сигнал u(t)', 'Исходный g(t)', "FontSize", 22);
xlim([-5, 5]);
ylim([-0.25 1.45]);
grid on;
%%
G = fftshift(fft(g));
U = fftshift(fft(u));
V = 1/dt;
dnu = 1/T;
nu = -V/2:dnu:V/2;
U_filt = U;
%%
n0 = 2;
U_filt(nu < -n0 | nu > n0 ) = 0;
%%
n1 = 0.5;
h = 0.1;
U_filt(nu < n1 + h & nu > n1 - h ) = 0;
U_filt(nu < -n1 + h & nu > -n1 - h ) = 0;
%%
hatnu = 5;
U_filt(nu > -hatnu & nu < hatnu) = 0;
%% Модули
f = figure();
plot(nu, abs(G), 'k', 'LineWidth', 12); hold on;
plot(nu, abs(U), 'r', 'LineWidth', 6);
plot(nu, abs(U_filt), 'g', 'LineWidth', 3);
legend('Модуль исходного', 'Модуль зашумлённого', 'Модуль отфильтрованного', 'FontSize', 21);
xlabel('Частота', "FontSize", 20); ylabel('Амплитуда', 'FontSize', 20);
title('Модули Фурье-образов', 'FontSize', 25);
xlim([-6-1 6+1]);
ylim([-10 650]);
grid on;
%%
u_filtered = ifft(ifftshift(U_filt));
%%
f = figure();
plot(t, u, 'b', 'LineWidth', 5); hold on;
plot(t, g, 'k--', 'LineWidth', 8);
plot(t, u_filtered, 'g', 'LineWidth', 4);
legend('Зашумленный сигнал u(t)', 'Исходный сигнал g(t)', 'Фильтрованный сигнал', "FontSize", 22);
xlabel('Время', 'FontSize', 20); ylabel('Значения', 'FontSize', 20);
title('Сравнение сигналов', 'FontSize', 25);
xlim([-5, 5]);
ylim([-0.6 1.45]);
grid on;
%%
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [100, 100, 1800, 1000]);
set(gca, 'FontSize', 22, 'GridLineWidth', 3, 'LineWidth', 2);
exportgraphics(f, "images\37.jpg")