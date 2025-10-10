% --- Настройки ---
clear; clc;

R = 2;
T = 1;
omega = @(n) 2*pi*n;

% Уровень аппроксимации
N = 10;
n_vals = -N:N;

% Сетка по времени
t = linspace(-1/8, 7/8, 2000);

% Определение f(t)
f = @(t) ...
    ((t >= -1/8 & t < 1/8) .* (2 + 16i * t) + ...
     (t >= 1/8  & t < 3/8) .* (4 - 16*t + 2i) + ...
     (t >= 3/8  & t < 5/8) .* (-2 + 1i*(8 - 16*t)) + ...
     (t >= 5/8  & t < 7/8) .* (-12 + 16*t - 2i)) ...
    .* (t >= -1/8 & t <= 7/8);

% Значения f(t) на сетке
ft = f(t);

% Для trapz-интегрирования
t_int = linspace(-1/8, 7/8, 2000);
ft_int = f(t_int);

% Вычисляем коэффициенты c_n с помощью trapz
c = zeros(size(n_vals));
for k = 1:length(n_vals)
    n = n_vals(k);
    integrand = ft_int .* exp(-1i * omega(n) * t_int);
    c(k) = trapz(t_int, integrand);
end

% Частичная сумма G_N(t)
GN = zeros(size(t));
for k = 1:length(n_vals)
    GN = GN + c(k) * exp(1i * omega(n_vals(k)) * t);
end

%% --- График Re(GN(t)) и Re(f(t)) ---
fi = figure();
plot(t, real(ft), 'k--', 'LineWidth', 5); hold on;
plot(t, real(GN), 'b', 'LineWidth', 6);
xlabel('t'); ylabel('Вещественная часть');
title(['Сравнение вещественных частей, N = ', num2str(N)]);
legend('Re(f(t))', ['Re(G_{', num2str(N), '}(t))'], "FontSize", 22, "Location", 'southwest');
grid on;

%% --- График Im(GN(t)) и Im(f(t)) ---
fi = figure();
plot(t, imag(ft), 'k--', 'LineWidth', 5); hold on;
plot(t, imag(GN), 'r', 'LineWidth', 6);
xlabel('t'); ylabel('Мнимая часть');
title(['Сравнение мнимых частей, N = ', num2str(N)]);
legend('Im(f(t))', ['Im(G_{', num2str(N), '}(t))'], "FontSize", 22);
grid on;
%%
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [100, 100, 1800, 1000]);
set(gca, 'FontSize', 22, 'GridLineWidth', 3, 'LineWidth', 2);
exportgraphics(fi, "images\37.jpg")
%%
% --- Проверка равенства Парсеваля для N = 10 ---
N = 10;
n_vals = -N:N;

% Сетка и значения f(t)
t_int = linspace(-1/8, 7/8, 2000);
ft_int = f(t_int);

% Левая часть: интеграл |f(t)|^2
LHS = trapz(t_int, abs(ft_int).^2);  % численная энергия

% Вычисление коэффициентов c_n
c = zeros(size(n_vals));
for k = 1:length(n_vals)
    n = n_vals(k);
    integrand = ft_int .* exp(-1i * omega(n) * t_int);
    c(k) = trapz(t_int, integrand);
end

% Правая часть: сумма квадратов модулей коэффициентов
RHS = sum(abs(c).^2);

% Вывод
fprintf('\n=== Проверка равенства Парсеваля при N = %d ===\n', N);
fprintf('Интеграл |f(t)|^2 (LHS)     = %.10f\n', LHS);
fprintf('Сумма |c_n|^2 (RHS, N=%d)   = %.10f\n', N, RHS);
fprintf('Абсолютная погрешность      = %.2e\n', abs(LHS - RHS));
