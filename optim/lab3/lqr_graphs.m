% График x(t) для всех систем на одном рисунке
figure;
plot(out.x, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояния объекта', 'Interpreter', 'tex');
legend('x_1', 'x_2', 'Location', 'best', 'Interpreter', 'tex');
saveas(gcf, 'images/x.png');

% График всех управлений u(t) на одном рисунке
figure;
plot(out.u, 'LineWidth', 2, 'DisplayName', 'u');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('u(t)', 'Interpreter', 'tex');
title('Оптимальное управление', 'Interpreter', 'tex');
legend('Location', 'best', 'Interpreter', 'tex');
saveas(gcf, 'images/u.png');


% График x(t) для всех систем на одном рисунке
figure;
plot(out.x1, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояния объекта', 'Interpreter', 'tex');
legend('x_1', 'x_2', 'Location', 'best', 'Interpreter', 'tex');
saveas(gcf, 'images/x1.png');

% График всех управлений u(t) на одном рисунке
figure;
plot(out.u1, 'LineWidth', 2, 'DisplayName', 'u');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('u(t)', 'Interpreter', 'tex');
title('Оптимальное управление', 'Interpreter', 'tex');
legend('Location', 'best', 'Interpreter', 'tex');
saveas(gcf, 'images/u1.png');


% График x(t) для всех систем на одном рисунке
figure;
plot(out.x2, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояния объекта', 'Interpreter', 'tex');
legend('x_1', 'x_2', 'Location', 'best', 'Interpreter', 'tex');
saveas(gcf, 'images/x2.png');

% График всех управлений u(t) на одном рисунке
figure;
plot(out.u2, 'LineWidth', 2, 'DisplayName', 'u');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('u(t)', 'Interpreter', 'tex');
title('Оптимальное управление', 'Interpreter', 'tex');
legend('Location', 'best', 'Interpreter', 'tex');
saveas(gcf, 'images/u2.png');