% --- A) Функция: вычисляет Q для заданной точки и всех сегментов ---
% Q_i = cos(psi_i)*(x - x_i) + sin(psi_i)*(y - y_i) - r_i
function Q = compute_Q_all(P, psi, r, x, y)
    % P: nx2, psi: (n-1)x1, r: (n-1)x1
    n = size(P,1);
    Q = zeros(n-1,1);
    for k = 1:n-1
        Q(k) = cos(psi(k))*(x - P(k,1)) + sin(psi(k))*(y - P(k,2)) - r(k);
    end
end

% --- B) Функция: возвращает активный сегмент для точки (x,y) ---
function idx = segment_id_by_Q(P, psi, r, x, y)
    Q = compute_Q_all(P, psi, r, x, y);
    % Правило: принадлежим первому i с Q(i) <= 0.
    idx = find(Q <= 0, 1, 'first');
    if isempty(idx)
        % Если не нашли (точка "впереди" всех сегментов) — возвращаем n-1
        idx = size(P,1)-1;
    end
end


% --- Визуализация полей Q_i ---
figure; hold on; axis equal; grid on;
plot(P(:,1), P(:,2), 'k--o', 'LineWidth', 1.5, 'MarkerFaceColor','k');
xlabel('x'); ylabel('y');
title('Области Q_i(x,y) \leq 0 для каждого сегмента');

% Сетка для оценки Q_i
[xg, yg] = meshgrid(linspace(min(P(:,1))-1, max(P(:,1))+1, 200), ...
                    linspace(min(P(:,2))-1, max(P(:,2))+1, 200));

% Для каждого сегмента построим границу Q_i=0
colors = lines(n-1);
for i = 1:n-1
    Qgrid = cos(psi(i))*(xg - P(i,1)) + sin(psi(i))*(yg - P(i,2)) - r(i);
    contour(xg, yg, Qgrid, [0 0], 'Color', colors(i,:), 'LineWidth', 1.2);
    contourf(xg, yg, Qgrid <= 0, 1, 'FaceColor', colors(i,:), ...
             'FaceAlpha', 0.1, 'EdgeColor', 'none');
    text(mean([P(i,1), P(i+1,1)]), mean([P(i,2), P(i+1,2)]), ...
         sprintf('Q_%d', i), 'FontSize', 9, 'FontWeight', 'bold', ...
         'Color', colors(i,:));
end

legend('Узлы', 'Location','bestoutside');