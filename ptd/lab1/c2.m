%% Нахождение траектории
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

% --- Задаём старт и финиш ---
start = sub2ind(size(map), 1, 1);   % верхний левый угол
goal  = sub2ind(size(map), 10, 10); % нижний правый угол

% --- Поиск пути с помощью A* ---
path = a_star(map, costs, start, goal);
if isempty(path)
    error('Путь не найден!');
end
mapSize = size(map);
[x, y] = ind2sub(mapSize, path);
x = flip(x);
y = flip(y);
P = [x', y'];

step = 0.0001;
r = 0.8;
k = 20;

pts_c2 = interpolate_c2(P, r, k, step);
pts_c1 = interpolate_c1(P, r, step);
pts_c0 = interpolate_c0(P, step);

%% Визуализация
a_star_plot(map, path, start, goal);
hold on;
plot(pts_c0(:,2), pts_c0(:,1), 'b-', 'LineWidth', 2);
plot(pts_c2(:,2), pts_c2(:,1), 'g-', 'LineWidth', 2);
hold off;
print('-djpeg', '-r600', 'images/c2.jpg');

%%
% pts_c0 = [x1 y1; x2 y2; ... ; xN yN];

x = pts_c2(:,1);
y = pts_c2(:,2);

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
xlim([0 1800]);
ylim([-0.05 1.4]);
%print('-djpeg', '-r600', 'images/c2_k.jpg');

%% Длина

% Вычисление расстояний между последовательными точками
distances_c2 = sqrt(diff(pts_c2(:,1)).^2 + diff(pts_c2(:,2)).^2);

% Суммируем все расстояния - это и есть длина кривой
curve_length_c2 = sum(distances_c2);

% Выводим результат
fprintf('Длина C2-гладкой кривой: %.3f единиц\n', curve_length_c2);

% Для наглядности можно построить график накопленной длины
cumulative_length_c2 = [0; cumsum(distances_c2)];

figure;
plot(cumulative_length_c2, 'LineWidth', 2, 'Color', [0.2, 0.6, 0.2]); % зеленый цвет
title('Накопленная длина C2-гладкой кривой');
xlabel('Номер точки (дискретизация)');
ylabel('Длина от начала');
xlim([0 1750]);
grid on;

% Показать финальную длину на графике
annotation('textbox', [0.6, 0.12, 0.1, 0.1], 'String', ...
    sprintf('Общая длина: %.3f', curve_length_c2), ...
    'BackgroundColor', 'white', 'FontSize', 10, 'EdgeColor', [0.2, 0.6, 0.2]);

print('-djpeg', '-r600', 'images/c2_length.jpg');