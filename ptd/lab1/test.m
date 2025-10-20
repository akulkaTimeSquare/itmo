% --- C1-сопряжённая траектория по формуле из изображения ---
clear; clc;

% --- входные данные ---
x = [0 2 4 6 8];    % координаты точек
y = [0 1 0 -1 0];
R = 1.0;            % радиус поворота

n = length(x);

% --- предварительные вычисления ---
psi = zeros(1, n-1);
r   = zeros(1, n-1);
for i = 1:n-1
    psi(i) = atan2(y(i+1) - y(i), x(i+1) - x(i));
    r(i) = sqrt((x(i+1) - x(i))^2 + (y(i+1) - y(i))^2);
end

% --- параметры дуг ---
sigma = zeros(1, n-2);
d_i   = zeros(1, n-2);
d_ci  = zeros(1, n-2);
delta = zeros(1, n-2);
xc    = zeros(1, n-2);
yc    = zeros(1, n-2);

for i = 1:n-2
    % угол между сегментами
    dpsi = psi(i+1) - psi(i);
    if dpsi > 0
        sigma(i) = pi - dpsi;
    else
        sigma(i) = -pi - dpsi;
    end

    % смещения и геометрия дуги
    d_i(i)  = abs(R / tan(sigma(i)/2));
    d_ci(i) = abs(R / sin(sigma(i)/2));
    delta(i) = pi - sigma(i)/2;

    % центр дуги
    Rlp = [cos(psi(i)) -sin(psi(i)); sin(psi(i)) cos(psi(i))];
    Rp  = [cos(delta(i)) -sin(delta(i)); sin(delta(i)) cos(delta(i))];
    center = [x(i+1); y(i+1)] + Rlp * Rp * [0; d_ci(i)];

    xc(i) = center(1);
    yc(i) = center(2);
end

% --- построение траектории через дуги ---
Sx = [];
Sy = [];

for i = 1:n-2
    theta_start = psi(i) - sigma(i)/2;
    theta_end   = psi(i) + sigma(i)/2;
    theta = linspace(theta_start, theta_end, 100);
    x_arc = xc(i) + R * cos(theta);
    y_arc = yc(i) + R * sin(theta);
    Sx = [Sx, x_arc];
    Sy = [Sy, y_arc];
end

% --- вывод ---
fprintf('C1 trajectory computed.\n');
disp('Arc centers [xc, yc]:');
disp([xc(:), yc(:)]);

% Результирующая траектория: Sx, Sy
% --- Визуализация C1-траектории ---
figure; hold on; axis equal; grid on;
title('C^1 Trajectory');
xlabel('x'); ylabel('y');

% исходные точки
plot(x, y, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Исходные точки');

% полученная траектория
plot(Sx, Sy, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Построенная траектория');

% центры дуг (если нужны для контроля)
if exist('xc','var')
    plot(xc, yc, 'ks', 'MarkerFaceColor','y', 'DisplayName', 'Центры дуг');
end

legend('Location','best');