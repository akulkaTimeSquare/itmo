function a_star_plot(map, costs, route, start, goal)
% A_STAR_PLOT — чистая визуализация карты и маршрута A*
% Без рамок, выровненная сетка, номера в центре клеток.

if ~islogical(map)
    error('map must be logical')
end

mapSize = size(map);

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

% --- Цвета ---
wallColor  = [0.25 0.25 0.25];
freeColor  = [0.92 0.92 0.92];
pathColor  = [0.65 0.3 0.7];
startColor = [0.35 0.8 0.35];
goalColor  = [0.3 0.45 0.9];

% --- normalize costs (если есть разные стоимости) ---
minCost = min(costs(:));
maxCost = max(costs(:));
mmCost = maxCost - minCost;
if mmCost == 0, mmCost = 1; end
costs = (costs(:) - minCost) / mmCost;

% --- Отрисовка карты ---
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

% --- Путь ---
[rp, cp] = ind2sub(mapSize, route);
plot(ha, cp, rp, ...
    'Color', pathColor, ...
    'LineWidth', 3, ...
    'LineStyle', '-');

% --- Старт и финиш ---
if nargin >= 4 && ~isempty(start)
    [rs, cs] = ind2sub(mapSize, start);
    rectangle('Position', [cs-0.5, rs-0.5, 1, 1], ...
        'FaceColor', startColor, ...
        'EdgeColor', 'none', ...
        'Curvature', 0.25);
end
if nargin >= 5 && ~isempty(goal)
    [rg, cg] = ind2sub(mapSize, goal);
    rectangle('Position', [cg-0.5, rg-0.5, 1, 1], ...
        'FaceColor', goalColor, ...
        'EdgeColor', 'none', ...
        'Curvature', 0.25);
end

% --- Настройки осей ---
xlim([0.5, mapSize(2)+0.5]);
ylim([0.5, mapSize(1)+0.5]);
xticks(1:mapSize(2));
yticks(1:mapSize(1));

ha.XTickLabel = string(1:mapSize(2));
ha.YTickLabel = string(1:mapSize(1));

hold(ha, 'off');
end