time = out.x.Time;
x = out.x.Data;
hatx = out.hatx.Data;
u = out.u.Data;
hatz = out.hatz.Data;

% Plot control input u
figure;
plot(time, u, 'LineWidth', 1.5);
xlabel('t');
ylabel('u(t)');
title('Управление системой для наблюдателя пониженного порядка');
grid on;
saveas(gcf, 'images/control_input.png');

% Plot z
figure;
plot(time, hatz, 'LineWidth', 1.5);
xlabel('t');
ylabel('ẑ');
title('Оценка скрытых состояний ẑ(t)');
if size(hatz,2) > 1
    legend_z = {};
    for i = 1:size(hatz,2)
        legend_z{end+1} = ['z_' num2str(i)];
    end
    legend(legend_z, 'Location', 'Best');
end
grid on;
saveas(gcf, 'images/z_plot.png');

% Plot states x and estimates hatx
n = size(x,2); % количество состояний
for i = 1:n
    figure;
    plot(time, x(:,i), 'b', 'LineWidth', 1.5); hold on;
    plot(time, hatx(:,i), 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('t');
    ylabel(['x_', num2str(i)]);
    % Легенда без использования \hat{}, чтобы избежать проблем отображения
    legend(['x_', num2str(i), '(t)'], ['x̂_', num2str(i), '(t)']);
    title(['Состояния x_', num2str(i), ' для наблюдателя пониженного порядка']);
    saveas(gcf, ['images/x', num2str(i), '_observer_down.png']);
end

% Plot estimation error
e = x - hatx;
figure;
hold on;
colors = lines(n);
for i = 1:n
    plot(time, e(:,i), 'Color', colors(i,:), 'LineWidth', 1.5);
end
grid on;
xlabel('t');
ylabel('e_i(t)');
legend(arrayfun(@(i) ['e_', num2str(i), '(t)'], 1:n, 'UniformOutput', false));
title('Ошибки оценок состояний e(t) = x(t) - x̂(t)');
saveas(gcf, 'images/estimation_error_observer_down.png');