T = 25;
dt = 0.001;
t = -T/2:dt:T/2;
p = zeros(size(t));
p(abs(t) <= 1/2) = 1;

nt = -T/2:0.001:T/2;
np = zeros(size(nt));
np(abs(nt) <= 1/2) = 1;

V = 100;
dnu = 0.01;
nu = -V/2:dnu:V/2;
nnu = -V/2:0.01:V/2;
hatp = sinc(nnu);

f = figure();

%plot(nt, p, 'b', 'LineWidth', 7);
%xlabel('Время', 'FontSize', 20); ylabel('Значения', "FontSize", 20);
%title('Прямоугольная функция П(t)', 'FontSize', 25);
%ylim([-0.1 1.1]);
%xlim([-1 1]);
%grid on;

%plot(nu, hatp, 'r', 'LineWidth', 8);
%xlabel('Частота', "FontSize", 20); ylabel('Амплитуда', 'FontSize', 20);
%title('Фурье-образ П(t) (аналитическое выражение)', 'FontSize', 25);
%ylim([-0.4 1.1]);
%xlim([-5 5]);
%grid on;


% Вычисление Фурье-образа численно
Pi_nu = zeros(size(nu));
for k = 1:length(nu)
    Pi_nu(k) = trapz(t, p .* exp(-1j*2*pi*nu(k)*t));
end

% Обратное преобразование Фурье
Pi_t_rec = zeros(size(t));
for k = 1:length(t)
    Pi_t_rec(k) = trapz(nu, Pi_nu .* exp(1j*2*pi*nu*t(k)));
end

% Нормировка (по аналогии с непрерывным интегралом)
Pi_t_rec = real(Pi_t_rec);

plot(t, Pi_t_rec, 'b', 'LineWidth', 7); hold on;
plot(nt, np, 'k--', LineWidth=6);
xlabel('Время', 'FontSize', 20); ylabel('Значения', "FontSize", 20);
title('Восстановленная прямоугольная волна (численное интегрирование)', 'FontSize', 25);
legend("Восстановленная функция", "Исходная функция", FontSize = 17);
ylim([-0.125 1.125]);
xlim([-1 1]);
grid on;

%plot(nu, Pi_nu, 'r', 'LineWidth', 9); hold on;
%plot(nnu, hatp, 'k', 'LineWidth', 4);
%xlabel('Частота', "FontSize", 20); ylabel('Амплитуда', 'FontSize', 20);
%title('Фурье-образ П(t) (численное интегрирование)', 'FontSize', 25);
%legend("Численный Фурье-образ", "Аналитический Фурье-образ", FontSize = 17)
%ylim([-0.4 1.1]);
%xlim([-V/2 V/2]);
%grid on;


set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [100, 100, 1800, 1000]);
set(gca, 'FontSize', 22, 'GridLineWidth', 3, 'LineWidth', 2);
%exportgraphics(f, "images\20.jpg")