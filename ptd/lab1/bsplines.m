clear; clc; close all;

% ===== Карта =====
map = logical([
    1 1 1 1 0 1 1 1 1 1;
    1 0 0 1 0 1 0 0 0 1;
    1 1 0 1 1 1 0 1 0 1;
    0 1 0 0 0 1 0 1 0 1;
    0 1 1 1 0 1 0 1 1 1;
    0 0 0 1 0 1 0 0 0 1;
    1 1 1 1 1 1 0 1 1 1;
    0 1 0 0 0 0 0 1 0 0;
    1 1 1 1 1 1 1 1 0 1;
    1 0 0 0 0 0 0 1 1 1;
]);

costs = ones(size(map));

% ===== Старт / финиш =====
start = sub2ind(size(map), 1, 1);
goal  = sub2ind(size(map), 10, 10);

% ===== Поиск пути A* =====
path = a_star(map, costs, start, goal);
if isempty(path)
    error('Путь не найден!');
end

mapSize = size(map);
[x, y] = ind2sub(mapSize, path);
x = flip(x);
y = flip(y);
P = [x', y'];

% ===== Параметры =====
step = 0.01;
degree = 3;

% ===== Построение B-spline =====
pts_c0 = interpolate_c0(P, step);
pts_bspline = interpolate_bsplines(P, degree, step);


%% Визуализация
a_star_plot(map, path, start, goal);
hold on;
plot(pts_c0(:,2), pts_c0(:,1), 'b-', 'LineWidth', 2);
plot(pts_bspline(:,2), pts_bspline(:,1), 'm-', 'LineWidth', 2);
hold off;
print('-djpeg', '-r600', 'images/bsplines.jpg');

%%
% pts_c0 = [x1 y1; x2 y2; ... ; xN yN];

x = pts_bspline(:,1);
y = pts_bspline(:,2);

N = length(x);

kappa = zeros(N,1);

% вычисляем локальные dt (расстояния между соседними точками)
dt_forward  = sqrt(diff(x).^2 + diff(y).^2);        % dt[i] = |P[i+1]-P[i]|
dt_backward = [dt_forward(1); dt_forward];         % для центральной разности

for i = 2:N-1
    dt_i = (dt_forward(i) + dt_backward(i))/2;  % усредненный шаг
    % первая производная
    dx = (x(i+1) - x(i-1)) / (2*dt_i);
    dy = (y(i+1) - y(i-1)) / (2*dt_i);
    % вторая производная
    ddx = (x(i+1) - 2*x(i) + x(i-1)) / (dt_i^2);
    ddy = (y(i+1) - 2*y(i) + y(i-1)) / (dt_i^2);
    % кривизна
    kappa(i) = abs(dx*ddy - dy*ddx) / (dx^2 + dy^2)^(3/2);
end

% концы
kappa(1) = kappa(2);
kappa(N) = kappa(N-1);


% визуализация
figure; 
plot(kappa, "LineWidth", 1.5);
title('Кривизна вдоль траектории'); 
xlabel('Точки'); 
ylabel('Кривизна'); 
grid on;
xlim([0 1300]);
ylim([-0.05 2.05])
print('-djpeg', '-r600', 'images/bsplines_k.jpg');