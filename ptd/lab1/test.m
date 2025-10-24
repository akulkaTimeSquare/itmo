%% Нахождение траектории
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
mapSize = size(map);
[x, y] = ind2sub(mapSize, path);

P = [x', y'];

% --- 4. Построение траектории ---
n = size(P,1);

R = 0.1;

% Предвычисляем psi_i и r_i для всех сегментов i=1..n-1
psi = zeros(n-1,1);
r   = zeros(n-1,1);
d   = zeros(n-1, 1);
d_c = zeros(n-1, 1);
delta = zeros(n-1, 1);
sigma = zeros(n-1, 1);
c = zeros(n-1, 2);
R_IP = @(ang) [cos(ang) sin(ang); -sin(ang) cos(ang)];

for i = 1:n-1
    dx = P(i+1,1) - P(i,1);
    dy = P(i+1,2) - P(i,2);
    psi(i) = atan2(dy, dx);
    r(i)   = sqrt(dx^2 + dy^2);
end

for i = 1:n-2
    if psi(i+1) - psi(i) > 0
        sigma(i) = pi - (psi(i+1) - psi(i));
    else
        sigma(i) = -pi - (psi(i+1) - psi(i));
    end
    d(i) = abs(R/tan(sigma(i)/2));
    d_c(i) = abs(R/sin(sigma(i)/2));
    c_i = P(i+1,:)' + R_IP(delta(i) + psi(i)) * [d_c(i); 0];
    c(i, 1) = c_i(1);
    c(i, 2) = c_i(2);
end

dt = 0.05;
v  = 0.25;
step = v*dt;
seg = 1;
pos = P(1,:);
traj = pos;

%%
max_steps = 10000;
mode = 'line';           % 'line' или 'arc'
for k = 1:max_steps
    
    % двигаемся вдоль направления psi(seg)
    pos = pos + step * [cos(psi(seg)), sin(psi(seg))]; 

    % проверяем Q для текущего сегмента:
    Qcur = cos(psi(seg))*(pos(1)-P(seg,1)) + sin(psi(seg))*(pos(2)-P(seg,2)) - r(seg);

    if Qcur > 0
        % пересекли границу сегмента -> переключаемся на следующий
        if seg < n-1
            seg = seg + 1;
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

%%
a_star_plot(map, path, start, goal);
hold on;
plot(c(:, 2), c(:, 1), 'o', ...
        'MarkerSize', 6, ...
        'LineWidth', 1.5);
hold off;


%% Визуализация
% --- Отрисовка траектории ---
a_star_plot(map, path, start, goal);
hold on;
plot(traj(:,2), traj(:,1), 'b-', 'LineWidth', 1.5);
hold off;
%print('-djpeg', '-r600', 'images/c1.jpg');