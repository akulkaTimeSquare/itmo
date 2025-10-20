clc; clear; close all

% --- 1. Фиксированная бинарная карта 10×10 ---
% 1 = проходимая клетка, 0 = стена
map = logical([
    1 1 1 1 0 1 1 1 1 1;
    1 0 0 1 0 1 0 0 0 1;
    1 1 0 1 1 1 0 1 0 1;
    0 1 0 0 0 1 0 1 0 1;
    0 1 1 1 0 1 0 1 1 1;
    0 0 0 1 0 1 0 0 0 1;
    1 1 1 1 1 1 0 1 1 1;
    1 0 0 0 0 0 0 1 0 0;
    1 1 1 1 1 1 1 1 0 1;
    1 0 0 0 0 0 0 1 1 1;
]);

costs = ones(size(map));

% --- 2. Задаём старт и финиш ---
start = sub2ind(size(map), 1, 1);   % верхний левый угол
goal  = sub2ind(size(map), 10, 10); % нижний правый угол

% --- 3. Поиск пути с помощью A* ---
path = a_star(map, costs, start, goal);

if isempty(path)
    error('Путь не найден!');
end

% --- 4. Подсчёт длины и числа поворотов ---
[r, c] = ind2sub(size(map), path);
len = numel(path);
dr = diff(r);
dc = diff(c);
dirs = atan2(dr, dc);
turns = sum(abs(diff(dirs)) > 1e-3);

fprintf('Длина пути: %d клеток, поворотов: %d\n', len, turns);

% --- 5. Визуализация ---
a_star_plot(map, costs, path);
title(sprintf('Фиксированная карта. Путь длиной %d и %d поворотов', len, turns));
print('-djpeg', '-r600', 'images/track.jpg');


%% c0
wallColor  = [0.25 0.25 0.25];
freeColor  = [0.92 0.92 0.92];
pathColor  = [0.65 0.3 0.7];
startColor = [0.35 0.8 0.35];
goalColor  = [0.9 0.3 0.3]; % красный для цели

hf = figure('Color', [1 1 1]);
ha = axes('Parent', hf, ...
    'YDir', 'reverse', ...
    'XDir', 'normal', ...
    'FontName', 'Segoe UI', ...
    'FontSize', 10, ...
    'Box', 'on', ...
    'Position', [0.02 0.08 0.96 0.92]);
grid(ha, 'on');
axis equal
hold(ha, 'on');

mapSize = size(map);
[x, y] = ind2sub(mapSize, path);

% --- Рисуем все клетки карты ---
for r = 1:mapSize(1)
    for c = 1:mapSize(2)
        if map(r, c)
            face = freeColor;
        else
            face = wallColor;
        end
        patch('Parent', ha, ...
            'XData', [c-0.5 c-0.5 c+0.5 c+0.5], ...
            'YData', [r-0.5 r+0.5 r+0.5 r-0.5], ...
            'FaceColor', face, ...
            'EdgeColor', [0.7 0.7 0.7], ...
            'LineWidth', 0.5);
    end
end

% --- Отрисовка траектории ---
n = length(x);
for i = 1:n-1
    % Параметрическая форма отрезка
    t = linspace(0, 1, 100);
    x_seg = x(i) + t * (x(i+1) - x(i));
    y_seg = y(i) + t * (y(i+1) - y(i));
    plot(y_seg, x_seg, 'b', 'LineWidth', 2);
end

% --- Кружок старта ---
[rs, cs] = ind2sub(mapSize, start);
plot(cs, rs, 'o', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', startColor, ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.5);

% --- Кружок цели ---
[rg, cg] = ind2sub(mapSize, goal);
plot(cg, rg, 'o', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', goalColor, ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.5);

xlim([0.5, mapSize(2)+0.5]);
ylim([0.5, mapSize(1)+0.5]);
xticks(1:mapSize(2));
yticks(1:mapSize(1));

ha.XTickLabel = string(1:mapSize(2));
ha.YTickLabel = string(1:mapSize(1));

hold(ha, 'off');