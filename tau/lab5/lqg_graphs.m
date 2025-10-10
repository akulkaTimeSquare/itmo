% График состояний и оценок для системы 1: (Q,R)
figure;
plot(out.x.Time, out.x.Data(:, 1), 'LineWidth', 2); 
hold on;
plot(out.xhat.Time, out.xhat.Data(:, 1), 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_1(t) и оценка x̂_1(t) при (Q_k,R_k) и (Q_l,R_l)', 'Interpreter', 'tex');
legend('x_1', 'x̂_1', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x_lqg_plot1.png');

figure;
plot(out.x.Time, out.x.Data(:, 2), 'LineWidth', 2); 
hold on;
plot(out.xhat.Time, out.xhat.Data(:, 2), 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_2(t) и оценка x̂_2(t) при (Q_k,R_k) и (Q_l,R_l)', 'Interpreter', 'tex');
legend('x_2', 'x̂_2', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x_lqg_plot2.png');

figure;
plot(out.x.Time, out.x.Data(:, 3), 'LineWidth', 2); 
hold on;
plot(out.xhat.Time, out.xhat.Data(:, 3), 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_3(t) и оценка x̂_3(t) при (Q_k,R_k) и (Q_l,R_l)', 'Interpreter', 'tex');
legend('x_3', 'x̂_3', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x_lqg_plot3.png');

figure;
plot(out.x.Time, out.x.Data(:, 4), 'LineWidth', 2); 
hold on;
plot(out.xhat.Time, out.xhat.Data(:, 4), 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_4(t) и оценка x̂_4(t) при (Q_k,R_k) и (Q_l,R_l)', 'Interpreter', 'tex');
legend('x_4', 'x̂_4', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x_lqg_plot4.png');

% График ошибки оценки для системы 1: (Q,R)
figure;
plot(out.xhat - out.x, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('e(t)', 'Interpreter', 'tex');
legend('e_1', 'e_2', 'e_3', 'e_4', 'Location', 'southeast', 'Interpreter', 'tex');
title('Ошибка наблюдателя e(t) при (Q_k,R_k) и (Q_l,R_l)', 'Interpreter', 'tex');
saveas(gcf, 'images/err_lqg_plot1.png');

% График всех управлений u(t) на одном рисунке
figure;
plot(out.u.Time, out.u.Data(:, 1), 'LineWidth', 2, 'DisplayName', 'u_1');
hold on;
plot(out.u.Time, out.u.Data(:, 2), 'LineWidth', 2, 'DisplayName', 'u_2');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('u(t)', 'Interpreter', 'tex');
title('Управления u(t)', 'Interpreter', 'tex');
legend('Location', 'best', 'Interpreter', 'tex');
saveas(gcf, 'images/controls_lqg.png');