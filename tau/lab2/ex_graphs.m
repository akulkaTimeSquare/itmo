time = out.x.Time;
x = out.x.Data;
hatx = out.hatx.Data;
u = out.u.Data;

% График управления
figure;
plot(time, u, 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('u(t)');
title('Управление u(t) при модальном управлении по выходу');
saveas(gcf, 'images/u_ex.png');


% Графики состояний x и оценок hatx
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
    title(['Состояние x_', num2str(i), ' и его оценка при модальном управлении по выходу']);
    saveas(gcf, ['images/x', num2str(i), '_ex.png']);
end

% График ошибок e = x - hatx (все на одном)
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
saveas(gcf, 'images/e_ex.png');