%% Пример: использование Q_i для выбора сегмента и симуляция движения
clear; clc; close all;

% Узловые точки (пример)
P = [0 0;
     2 1;
     4 1.5;
     6 0];

n = size(P,1);

% Предвычисляем psi_i и r_i для всех сегментов i=1..n-1
psi = zeros(n-1,1);
r   = zeros(n-1,1);
for i = 1:n-1
    dx = P(i+1,1) - P(i,1);
    dy = P(i+1,2) - P(i,2);
    psi(i) = atan2(dy, dx);
    r(i)   = sqrt(dx^2 + dy^2);
end

dt = 0.05;
v  = 0.25;
step = v*dt;
seg = 1;
pos = P(1,:);
traj = pos;

max_steps = 10000;
for k = 1:max_steps
    % двигаемся вдоль направления psi(seg)
    pos = pos + step * [cos(psi(seg)), sin(psi(seg))]; 

    % проверяем Q для текущего сегмента:
    Qcur = cos(psi(seg))*(pos(1)-P(seg,1)) + sin(psi(seg))*(pos(2)-P(seg,2)) - r(seg);

    if Qcur > 0
        % пересекли границу сегмента -> переключаемся на следующий
        if seg < n-1
            seg = seg + 1;
            % Опционально: корректируем pos, чтобы не "залипать" за границей
            % (сдвигаем немного назад вдоль нового сегмента)
            % pos = P(seg,:) + 0.001 * [cos(psi(seg)), sin(psi(seg))];
        else
            % дошли до последнего сегмента и вышли за него
            traj = [traj; pos];
            break;
        end
    end

    traj = [traj; pos];

    % стоп — если прошли последний узел (по расстоянию)
    if norm(pos - P(end,:)) < 1e-3
        break;
    end
end

% --- Визуализация ---
figure; hold on; axis equal; grid on;
plot(P(:,1), P(:,2), 'k--o', 'LineWidth', 1.2, 'MarkerFaceColor','k');
plot(traj(:,1), traj(:,2), 'b-', 'LineWidth', 1.5);
xlabel('x'); ylabel('y');
title('Симуляция движения с переключением сегмента по Q_i');
legend('Узлы','Траектория агента','Тестовые точки','Location','best');
