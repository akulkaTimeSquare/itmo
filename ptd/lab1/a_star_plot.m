function a_star_plot(map, route, start, goal)
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
    goalColor  = [0.8 0.3 0.3];

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
    rp = rp(2:length(rp)-1);
    cp = cp(2:length(cp)-1);
    plot(ha, cp, rp, 'o', ...
        'MarkerFaceColor', pathColor, ...
        'MarkerEdgeColor', pathColor, ...
        'MarkerSize', 6, ...
        'LineWidth', 1.5);

    % --- Старт и финиш ---
    if nargin >= 4 && ~isempty(start)
        [rs, cs] = ind2sub(mapSize, start);
        plot(cs, rs, 'o', ...
            'MarkerSize', 12, ...
            'MarkerFaceColor', startColor, ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.5);
    end

    if nargin >= 4 && ~isempty(goal)
        [rg, cg] = ind2sub(mapSize, goal);
        plot(cg, rg, 'o', ...
            'MarkerSize', 12, ...
            'MarkerFaceColor', goalColor, ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.5);
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