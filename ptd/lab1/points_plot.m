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

% --- 2. Задаём старт и финиш ---
start = sub2ind(size(map), 1, 1);   % верхний левый угол
goal  = sub2ind(size(map), 10, 10); % нижний правый угол

% --- 3. Поиск пути с помощью A* ---
path = a_star(map, costs, start, goal);

if isempty(path)
    error('Путь не найден!');
end


% --- 5. Визуализация ---
a_star_plot(map, path, start, goal);
print('-djpeg', '-r600', 'images/points_plot.jpg');