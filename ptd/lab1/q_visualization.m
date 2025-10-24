% --- A) Функция: вычисляет Q для заданной точки и всех сегментов ---
function Q = compute_Q_all(P, psi, r, x, y)
    n = size(P,1);
    Q = zeros(n-1,1);
    for k = 1:n-1
        Q(k) = cos(psi(k))*(x - P(k,1)) + sin(psi(k))*(y - P(k,2)) - r(k);
    end
end

% --- B) Функция: возвращает активный сегмент для точки (x,y) ---
function idx = segment_id_by_Q(P, psi, r, x, y)
    Q = compute_Q_all(P, psi, r, x, y);
    idx = find(Q <= 0, 1, 'first');
    if isempty(idx)
        idx = size(P,1)-1;
    end
end

% --- Исходные данные ---
P = [0 0;
     2 0;
     2 2];          % точки траектории
n = size(P,1);

R = 0.5;            % радиус сглаживания дуг
v = 0.25;           % скорость
dt = 0.05;          % шаг по времени
step = v*dt;

% --- Предвычисления углов и длин ---
psi = zeros(n-1,1);
r   = zeros(n-1,1);
for i = 1:n-1
    d = P(i+1,:) - P(i,:);
    psi(i) = atan2(d(2), d(1));
    r(i) = norm(d);
end

% --- УЛУЧШЕННАЯ ВИЗУАЛИЗАЦИЯ ---
figure('Position', [100, 100, 1200, 900]);
hold on; axis equal; grid on;

% Цветовая схема
colors = lines(n-1);
dark_colors = colors * 0.7; % Более темные версии для контуров

% Сетка для оценки Q_i с более высоким разрешением
[xg, yg] = meshgrid(linspace(min(P(:,1))-1.5, max(P(:,1))+1.5, 300), ...
                    linspace(min(P(:,2))-1.5, max(P(:,2))+1.5, 300));

% Создаем маску для каждой области Q_i <= 0
region_masks = false(size(xg, 1), size(xg, 2), n-1);

% Вычисляем области и строим красивые контуры
for i = 1:n-1
    Qgrid = cos(psi(i))*(xg - P(i,1)) + sin(psi(i))*(yg - P(i,2)) - r(i);
    region_masks(:,:,i) = Qgrid <= 0;
    
    % Красивый контур с плавными линиями
    contour(xg, yg, Qgrid, [0 0], ...
           'Color', dark_colors(i,:), ...
           'LineWidth', 2.5, ...
           'LineStyle', '-');
    
    % Заливка областей с прозрачностью
    contourf(xg, yg, region_masks(:,:,i), 1, ...
            'FaceColor', colors(i,:), ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
end

% Визуализация траектории
plot(P(:,1), P(:,2), 'k-', 'LineWidth', 3, 'Color', [0.2 0.2 0.2]);
scatter(P(:,1), P(:,2), 120, 'filled', ...
        'MarkerFaceColor', [0.9 0.1 0.1], ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 2);

% Подписи узлов
for i = 1:size(P,1)
    text(P(i,1)+0.1, P(i,2)+0.1, sprintf('P%d', i), ...
         'FontSize', 12, 'FontWeight', 'bold', ...
         'Color', [0.7 0 0], 'BackgroundColor', 'white');
end

% Стрелки направления для сегментов
for i = 1:n-1
    mid_point = (P(i,:) + P(i+1,:)) / 2;
    arrow_length = 0.3;
    
    % Вычисляем направление стрелки
    dx = cos(psi(i)) * arrow_length;
    dy = sin(psi(i)) * arrow_length;
    
    quiver(mid_point(1), mid_point(1), dx, dy, ...
           'Color', colors(i,:), ...
           'LineWidth', 2, ...
           'MaxHeadSize', 1.5, ...
           'AutoScale', 'off');
end

% Создаем легенду
legend_entries = cell(n,1);
legend_entries{1} = 'Траектория';
for i = 1:n-1
    legend_entries{i+1} = sprintf('Область Q_%d ≤ 0', i);
end

% Отображение областей сегментов в легенде
for i = 1:n-1
    plot(NaN, NaN, '-', 'Color', colors(i,:), 'LineWidth', 3);
end

legend(legend_entries, 'Location', 'northeastoutside', ...
       'FontSize', 10, 'Box', 'off');

% Настройка осей и заголовка
xlabel('Координата X', 'FontSize', 12);
ylabel('Координата Y', 'FontSize', 12);
title('Визуализация областей Q_i(x,y) ≤ 0 для сегментов траектории', ...
      'FontSize', 14, 'FontWeight', 'bold');

% Добавляем информационную панель
info_str = {sprintf('Точки траектории: %d', n), ...
           sprintf('Сегменты: %d', n-1), ...
           sprintf('Радиус сглаживания: %.1f', R), ...
           sprintf('Скорость: %.2f', v)};

annotation('textbox', [0.02, 0.02, 0.3, 0.15], ...
           'String', info_str, ...
           'FontSize', 10, ...
           'BackgroundColor', [0.95 0.95 0.95], ...
           'EdgeColor', [0.5 0.5 0.5]);

% Улучшаем внешний вид графика
set(gca, 'FontSize', 11, ...
         'GridAlpha', 0.3, ...
         'Box', 'on');

% Добавляем цветовую карту для отображения активных сегментов
figure('Position', [100, 100, 1000, 800]);
hold on; axis equal; grid on;

% Вычисляем активные сегменты для каждой точки сетки
active_segments = zeros(size(xg));
for i = 1:numel(xg)
    active_segments(i) = segment_id_by_Q(P, psi, r, xg(i), yg(i));
end

% Визуализация карты активных сегментов
imagesc(xg(1,:), yg(:,1), active_segments);
colormap(colors);
colorbar('Ticks', 1:n-1, 'TickLabels', arrayfun(@(x) sprintf('Seg %d', x), 1:n-1, 'UniformOutput', false));

% Наложение траектории
plot(P(:,1), P(:,2), 'k-', 'LineWidth', 4);
scatter(P(:,1), P(:,2), 100, 'filled', ...
        'MarkerFaceColor', 'white', ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 2);
xlim([-1.5 3.5]);

xlabel('Координата X', 'FontSize', 12);
ylabel('Координата Y', 'FontSize', 12);
title('Карта активных сегментов по правилу Q_i ≤ 0', ...
      'FontSize', 14, 'FontWeight', 'bold');

set(gca, 'FontSize', 11, 'YDir', 'normal');