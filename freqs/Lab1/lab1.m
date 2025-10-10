% Параметры
T = 2;                         % Период
omega0 = 2*pi/T;               % Основная частота
N = 100;                        % Число гармоник
Nt = 1000;                     % Точек для интегрирования
tb = linspace(1, 3, Nt);       % Один период: [1, 3)
dt = tb(2) - tb(1);

% Определение f(t) на одном периоде
f_base = ones(size(tb));       % [1, 2)
f_base(tb >= 2) = 2;           % [2, 3)

% Расчет коэффициентов на одном периоде
a0 = (2/T) * trapz(tb, f_base);    % Средняя составляющая

a = zeros(1, N);
b = zeros(1, N);
for n = 1:N
    a(n) = (2/T) * trapz(tb, f_base .* cos(n * omega0 * tb));
    b(n) = (2/T) * trapz(tb, f_base .* sin(n * omega0 * tb));
end

n_range = -N:N;                               % Индексы для c_n
c = zeros(size(n_range));                    % Массив комплексных коэффициентов
for k = 1:length(n_range)
    n = n_range(k);
    c(k) = (1/T) * trapz(tb, f_base .* exp(-1i * n * omega0 * tb));
end


t_start = 1;
t_end = t_start+3*T+1;
% Восстановление функции на более широком интервале (3 периода)
t = linspace(t_start, t_end, 2000);          % Интервал от 1 до 7
f_approx = a0/2 * ones(size(t));  % Постоянная часть

for n = 1:N
    f_approx = f_approx + a(n)*cos(n * omega0 * t) + b(n)*sin(n * omega0 * t);
end

f_complex = zeros(size(t));
for k = 1:length(n_range)
    n = n_range(k);
    f_complex = f_complex + c(k) * exp(1i * n * omega0 * t);
end
f_complex = real(f_complex);

% Оригинальная функция на интервале [1, 7]
t_mod = mod(t - 1, T);
f_true = ones(size(t));
f_true(t_mod >= 1) = 2;

E_func = (1/T) * trapz(tb, f_base.^2)
E_real = (a0^2)/4 + 0.5*sum(a.^2 + b.^2)
E_complex = sum(abs(c).^2)



f = figure();
plot(t, f_true, 'k--', 'LineWidth', 5); hold on;
plot(t, f_approx, 'b', 'LineWidth', 6);
plot(t, f_complex, 'r--', 'LineWidth', 6);
xlabel('Аргумент', 'FontSize', 20); ylabel('Значения', "FontSize", 20);
%title(['Квадратная волна'], 'FontSize', 25);
title(['Ряды Фурье при N = ' num2str(N)], 'FontSize', 25);
legend('Исходная функция', 'Приближение вещественным рядом', ...
    'Приближение комплексным рядом', "FontSize", 20);
ylim([0.5, 2.5])
grid on;

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [100, 100, 1800, 1000]);
set(gca, 'FontSize', 22, 'GridLineWidth', 3, 'LineWidth', 2);
%exportgraphics(f, "images\8.jpg");