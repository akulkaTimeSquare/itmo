function a_star_plot(map, route, start, goal)
    % A_STAR_PLOT — визуализация карты и маршрута A* (версия для отчёта, улучшенный контраст)

    if ~islogical(map)
        error('map must be logical')
    end

    mapSize = size(map);

    hf = figure('Color', [1 1 1]); % белый фон
    ha = axes('Parent', hf, ...
        'YDir', 'reverse', ...
        'XDir', 'normal', ...
        'FontName', 'Segoe UI', ...
        'FontSize', 10, ...
        'Box', 'on', ...
        'Color', [1 1 1], ...
        'XColor', [0 0 0], ...
        'YColor', [0 0 0], ...
        'Position', [0.05 0.08 0.93 0.9]);
    grid(ha, 'on');
    axis equal
    hold(ha, 'on');

    % --- Цветовая схема (отчёт + чуть больше контраста) ---
    wallColor  = [0.55 0.55 0.55];   % темно-серый — препятствия (раньше было 0.75)
    freeColor  = [1.0 1.0 1.0];      % белый — свободные клетки
    pathColor  = [0.25 0.35 0.85];   % глубокий синий — путь
    startColor = [0.2 0.65 0.3];     % зелёный — старт
    goalColor  = [0.85 0.25 0.25];   % красный — цель
    gridColor  = [0.85 0.85 0.85];   % чуть светлее, чтобы не мешала серым стенам

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
                'EdgeColor', gridColor, ...
                'LineWidth', 0.4, ...
                'HandleVisibility', 'off');
        end
    end

    % --- Путь ---
    [rp, cp] = ind2sub(mapSize, route);
    rp = rp(1:end);
    cp = cp(1:end);
    plot(ha, cp, rp, 'o', ...
        'Color', pathColor, ...
        'LineWidth', 1.8, ...
        'HandleVisibility', 'off');

    % --- Старт ---
    if nargin >= 3 && ~isempty(start)
        [rs, cs] = ind2sub(mapSize, start);
        plot(cs, rs, 'o', ...
            'MarkerSize', 10, ...
            'MarkerFaceColor', startColor, ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 0.8, ...
            'HandleVisibility', 'off');
    end

    % --- Цель ---
    if nargin >= 4 && ~isempty(goal)
        [rg, cg] = ind2sub(mapSize, goal);
        plot(cg, rg, 'o', ...
            'MarkerSize', 10, ...
            'MarkerFaceColor', goalColor, ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 0.8, ...
            'HandleVisibility', 'off');
    end

    % --- Настройки осей ---
    xlim([0.5, mapSize(2)+0.5]);
    ylim([0.5, mapSize(1)+0.5]);
    xticks(1:mapSize(2));
    yticks(1:mapSize(1));
    ha.XTickLabel = string(1:mapSize(2));
    ha.YTickLabel = string(1:mapSize(1));

    grid(ha, 'on');
    ha.GridColor = gridColor;
    ha.GridAlpha = 0.5;

    hold(ha, 'off');
end
