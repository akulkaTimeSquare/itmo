%% 10
figure;
plot(out.x, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояния x(t) при \gamma = 10', 'Interpreter', 'tex');
legend('x_1', 'x_2', 'Interpreter', 'tex');
ylim([-0.35 0.2]);
saveas(gcf, 'images/x_10.png');

figure;
plot(out.u, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('u(t)', 'Interpreter', 'tex');
title('График управления u(t) при \gamma = 10', 'Interpreter', 'tex');
legend('u', 'Interpreter', 'tex');
saveas(gcf, 'images/u_10.png');

%% 100
figure;
plot(out.x1, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояния x(t) при \gamma = 100', 'Interpreter', 'tex');
legend('x_1', 'x_2', 'Interpreter', 'tex');
ylim([-0.21, 0.16]);
saveas(gcf, 'images/x_100.png');

figure;
plot(out.u1, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('u(t)', 'Interpreter', 'tex');
title('График управления u(t) при \gamma = 100', 'Interpreter', 'tex');
legend('u', 'Interpreter', 'tex');
saveas(gcf, 'images/u_100.png');

%% diff
figure;
plot(out.x1, 'LineWidth', 2); hold on;
plot(out.x, 'LineWidth', 2, 'LineStyle', '--'); 
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Сравнение состояний x(t) при разных \gamma', 'Interpreter', 'tex');
legend('x_1 при \gamma = 100', 'x_2 при \gamma = 100', 'x_1 при \gamma = 10', 'x_2 при \gamma = 10', 'Interpreter', 'tex');
ylim([-0.35 0.2]);
saveas(gcf, 'images/x_diff.png');

figure;
plot(out.u1, 'LineWidth', 2); hold on;
plot(out.u, 'LineWidth', 2, 'LineStyle', '--'); 
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('u(t)', 'Interpreter', 'tex');
title('Сравнение управлений u(t) при разных \gamma', 'Interpreter', 'tex');
legend('u при \gamma = 100', 'u при \gamma = 10', 'Interpreter', 'tex');
saveas(gcf, 'images/u_diff.png');