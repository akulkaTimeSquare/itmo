% График x1(t) для всех систем на одном рисунке
figure;
plot(out.x1, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояния x(t) при (Q,R)', 'Interpreter', 'tex');
legend('x_1', 'x_2', 'x_3', 'Location', 'best', 'Interpreter', 'tex');
saveas(gcf, 'images/x1_lqr_plot1.png');

% График x2(t) для всех систем на отдельном рисунке
figure;
plot(out.x2, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояния x(t) при (aQ,R)', 'Interpreter', 'tex');
legend('x_1', 'x_2', 'x_3', 'Location', 'best', 'Interpreter', 'tex');
saveas(gcf, 'images/x2_lqr_plot2.png');

% График x3(t) для всех систем на отдельном рисунке
figure;
plot(out.x3, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояния x(t) при (Q,aR)', 'Interpreter', 'tex');
legend('x_1', 'x_2', 'x_3', 'Location', 'best', 'Interpreter', 'tex');
saveas(gcf, 'images/x3_lqr_plot3.png');

% График x4(t) для всех систем на отдельном рисунке
figure;
plot(out.x4, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояния x(t) при (aQ,aR)', 'Interpreter', 'tex');
legend('x_1', 'x_2', 'x_3', 'Location', 'best', 'Interpreter', 'tex');
saveas(gcf, 'images/x4_lqr_plot4.png');

% График всех управлений u(t) на одном рисунке
figure;
plot(out.u1, 'LineWidth', 2, 'DisplayName', '(Q,R)');
hold on;
plot(out.u2, 'LineWidth', 2, 'DisplayName', '(aQ,R)');
plot(out.u3, 'LineWidth', 2, 'DisplayName', '(Q,aR)');
plot(out.u4, 'LineWidth', 2, 'DisplayName', '(aQ,aR)', 'LineStyle', '--', 'Color', "k");
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('u(t)', 'Interpreter', 'tex');
title('Управления u(t)', 'Interpreter', 'tex');
legend('Location', 'best', 'Interpreter', 'tex');
saveas(gcf, 'images/controls_lqr.png');

