% Fourier coefficients using trapezoidal rule + partial sums visualization
clear; clc;

% Parameters
R = 2;
T = 1;
omega = @(n) 2*pi*n;

% Time discretization for integration
Nt = 100000;
t_int = linspace(-1/8, 7/8, Nt);
dt = t_int(2) - t_int(1);

% Piecewise function f(t), zero outside domain
func = @(t) ...
    ((t >= -1/8 & t < 1/8) .* (2 + 16i * t) + ...
     (t >= 1/8  & t < 3/8) .* (4 - 16*t + 2i) + ...
     (t >= 3/8  & t < 5/8) .* (-2 + 1i*(8 - 16*t)) + ...
     (t >= 5/8  & t < 7/8) .* (-12 + 16*t - 2i)) ...
    .* (t >= -1/8 & t <= 7/8);

% Precompute f(t) on integration grid
ft = func(t_int);

% Time vector for plotting (including closure point)
t = linspace(-1/8, 7/8, 2000);
t = [t, t(1)];
ft_plot = func(t);
valid = (ft_plot ~= 0);  % mask to remove spurious points

% Plot f(t)
f = figure();
plot(real(ft_plot(valid)), imag(ft_plot(valid)), 'k', 'LineWidth', 6); hold on;

% Evaluate and plot G_N(t)
Ns = [1, 2, 3, 10];
colors = ['r', 'g', 'b', 'm'];
widths = [12, 6, 7, 8];

for i = 1:length(Ns)
    N = Ns(i);
    n_vals = -N:N;
    c = zeros(size(n_vals));
    
    % Compute Fourier coefficients using trapezoidal rule
    for k = 1:length(n_vals)
        n = n_vals(k);
        c(k) = trapz(t_int, ft .* exp(-1i * omega(n) * t_int));  % numerical integration over t_int
    end

    % Compute partial sum G_N(t)
    GN = zeros(size(t));
    for k = 1:length(n_vals)
        GN = GN + c(k) * exp(1i * omega(n_vals(k)) * t);
    end

    % Plot partial sum
    plot(real(GN), imag(GN), colors(i), 'LineWidth', widths(i));
    
    % Report coefficients for N = 2
    if N == 2
        disp('--- Коэффициенты c_n при N = 2 (trapz) ---');
        disp(table(n_vals.', real(c(:)), imag(c(:)), abs(c(:)), ...
            'VariableNames', {'n', 'Re_cn', 'Im_cn', 'Abs_cn'}));
    end
end

legend('f(t)', 'G_1(t)', 'G_2(t)', 'G_3(t)', 'G_{10}(t)', "FontSize", 22);
xlabel('Вещественная часть', 'FontSize', 20); ylabel('Мнимая часть', 'FontSize', 20);
xlim([-2.5 2.5]);
ylim([-2.5, 2.5]);
title('Параметрические графики f(t) и G_N(t) на комплексной плоскости', 'FontSize', 25);
axis equal; grid on;
%%
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [100, 100, 1800, 1000]);
set(gca, 'FontSize', 22, 'GridLineWidth', 3, 'LineWidth', 2);
exportgraphics(f, "images\30.jpg")
