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
a_star_plot(map, path, start, goal);
title(sprintf('Фиксированная карта. Путь длиной %d и %d поворотов', len, turns));
print('-djpeg', '-r600', 'images/points_plot.jpg');