% График состояний и оценок для системы 1: (Q,R)
figure;
plot(out.x1.Time, out.x1.Data(:, 1), 'LineWidth', 2); 
hold on;
plot(out.xhat1.Time, out.xhat1.Data(:, 1), 'LineWidth', 2, 'LineStyle', '-.');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_1(t) и оценка x̂_1(t) при (Q,R)', 'Interpreter', 'tex');
legend('x_1', 'x̂_1', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x1_1_lqe_plot1.png');

figure;
plot(out.x1.Time, out.x1.Data(:, 2), 'LineWidth', 2); 
hold on;
plot(out.xhat1.Time, out.xhat1.Data(:, 2), 'LineWidth', 2, 'LineStyle', '-.');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_2(t) и оценка x̂_2(t) при (Q,R)', 'Interpreter', 'tex');
legend('x_2', 'x̂_2', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x1_2_lqe_plot1.png');

figure;
plot(out.x1.Time, out.x1.Data(:, 3), 'LineWidth', 2); 
hold on;
plot(out.xhat1.Time, out.xhat1.Data(:, 3), 'LineWidth', 2, 'LineStyle', '-.');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_3(t) и оценка x̂_3(t) при (Q,R)', 'Interpreter', 'tex');
legend('x_3', 'x̂_3', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x1_3_lqe_plot1.png');

figure;
plot(out.x1.Time, out.x1.Data(:, 4), 'LineWidth', 2); 
hold on;
plot(out.xhat1.Time, out.xhat1.Data(:, 4), 'LineWidth', 2, 'LineStyle', '-.');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_4(t) и оценка x̂_4(t) при (Q,R)', 'Interpreter', 'tex');
legend('x_4', 'x̂_4', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x1_4_lqe_plot1.png');

% График ошибки оценки для системы 1: (Q,R)
figure;
plot(out.xhat1 - out.x1, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('e(t)', 'Interpreter', 'tex');
legend('e_1', 'e_2', 'e_3', 'e_4', 'Location', 'southeast', 'Interpreter', 'tex');
title('Ошибка наблюдателя e(t) при (Q,R)', 'Interpreter', 'tex');
saveas(gcf, 'images/err1_lqe_plot1.png');

% График состояний и оценок для системы 2: (aQ,R)
figure;
plot(out.x2.Time, out.x2.Data(:, 1), 'LineWidth', 2); 
hold on;
plot(out.xhat2.Time, out.xhat2.Data(:, 1), 'LineWidth', 2, 'LineStyle', '-.');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_1(t) и оценка x̂_1(t) при (aQ,R)', 'Interpreter', 'tex');
legend('x_1', 'x̂_1', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x2_1_lqe_plot2.png');

figure;
plot(out.x2.Time, out.x2.Data(:, 2), 'LineWidth', 2); 
hold on;
plot(out.xhat2.Time, out.xhat2.Data(:, 2), 'LineWidth', 2, 'LineStyle', '-.');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_2(t) и оценка x̂_2(t) при (aQ,R)', 'Interpreter', 'tex');
legend('x_2', 'x̂_2', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x2_2_lqe_plot2.png');

figure;
plot(out.x2.Time, out.x2.Data(:, 3), 'LineWidth', 2); 
hold on;
plot(out.xhat2.Time, out.xhat2.Data(:, 3), 'LineWidth', 2, 'LineStyle', '-.');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_3(t) и оценка x̂_3(t) при (aQ,R)', 'Interpreter', 'tex');
legend('x_3', 'x̂_3', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x2_3_lqe_plot2.png');

figure;
plot(out.x2.Time, out.x2.Data(:, 4), 'LineWidth', 2); 
hold on;
plot(out.xhat2.Time, out.xhat2.Data(:, 4), 'LineWidth', 2, 'LineStyle', '-.');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_4(t) и оценка x̂_4(t) при (aQ,R)', 'Interpreter', 'tex');
legend('x_4', 'x̂_4', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x2_4_lqe_plot2.png');

% График ошибки оценки для системы 2: (aQ,R)
figure;
plot(out.xhat2 - out.x2, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('e(t)', 'Interpreter', 'tex');
legend('e_1', 'e_2', 'e_3', 'e_4', 'Location', 'southeast', 'Interpreter', 'tex');
title('Ошибка наблюдателя e(t) при (aQ,R)', 'Interpreter', 'tex');
saveas(gcf, 'images/err2_lqe_plot2.png');

% График состояний и оценок для системы 3: (Q,aR)
figure;
plot(out.x3.Time, out.x3.Data(:, 1), 'LineWidth', 2); 
hold on;
plot(out.xhat3.Time, out.xhat3.Data(:, 1), 'LineWidth', 2, 'LineStyle', '-.');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_1(t) и оценка x̂_1(t) при (Q,aR)', 'Interpreter', 'tex');
legend('x_1', 'x̂_1', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x3_1_lqe_plot3.png');

figure;
plot(out.x3.Time, out.x3.Data(:, 2), 'LineWidth', 2); 
hold on;
plot(out.xhat3.Time, out.xhat3.Data(:, 2), 'LineWidth', 2, 'LineStyle', '-.');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_2(t) и оценка x̂_2(t) при (Q,aR)', 'Interpreter', 'tex');
legend('x_2', 'x̂_2', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x3_2_lqe_plot3.png');

figure;
plot(out.x3.Time, out.x3.Data(:, 3), 'LineWidth', 2); 
hold on;
plot(out.xhat3.Time, out.xhat3.Data(:, 3), 'LineWidth', 2, 'LineStyle', '-.');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_3(t) и оценка x̂_3(t) при (Q,aR)', 'Interpreter', 'tex');
legend('x_3', 'x̂_3', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x3_3_lqe_plot3.png');

figure;
plot(out.x3.Time, out.x3.Data(:, 4), 'LineWidth', 2); 
hold on;
plot(out.xhat3.Time, out.xhat3.Data(:, 4), 'LineWidth', 2, 'LineStyle', '-.');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_4(t) и оценка x̂_4(t) при (Q,aR)', 'Interpreter', 'tex');
legend('x_4', 'x̂_4', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x3_4_lqe_plot3.png');

% График ошибки оценки для системы 3: (Q,aR)
figure;
plot(out.xhat3 - out.x3, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('e(t)', 'Interpreter', 'tex');
legend('e_1', 'e_2', 'e_3', 'e_4', 'Location', 'southeast', 'Interpreter', 'tex');
title('Ошибка наблюдателя e(t) при (Q,aR)', 'Interpreter', 'tex');
saveas(gcf, 'images/err3_lqe_plot3.png');

% График состояний и оценок для системы 4: (aQ,aR)
figure;
plot(out.x4.Time, out.x4.Data(:, 1), 'LineWidth', 2); 
hold on;
plot(out.xhat4.Time, out.xhat4.Data(:, 1), 'LineWidth', 2, 'LineStyle', '-.');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_1(t) и оценка x̂_1(t) при (aQ,aR)', 'Interpreter', 'tex');
legend('x_1', 'x̂_1', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x4_1_lqe_plot4.png');

figure;
plot(out.x4.Time, out.x4.Data(:, 2), 'LineWidth', 2); 
hold on;
plot(out.xhat4.Time, out.xhat4.Data(:, 2), 'LineWidth', 2, 'LineStyle', '-.');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_2(t) и оценка x̂_2(t) при (aQ,aR)', 'Interpreter', 'tex');
legend('x_2', 'x̂_2', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x4_2_lqe_plot4.png');

figure;
plot(out.x4.Time, out.x4.Data(:, 3), 'LineWidth', 2); 
hold on;
plot(out.xhat4.Time, out.xhat4.Data(:, 3), 'LineWidth', 2, 'LineStyle', '-.');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_3(t) и оценка x̂_3(t) при (aQ,aR)', 'Interpreter', 'tex');
legend('x_3', 'x̂_3', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x4_3_lqe_plot4.png');

figure;
plot(out.x4.Time, out.x4.Data(:, 4), 'LineWidth', 2); 
hold on;
plot(out.xhat4.Time, out.xhat4.Data(:, 4), 'LineWidth', 2, 'LineStyle', '-.');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояние x_4(t) и оценка x̂_4(t) при (aQ,aR)', 'Interpreter', 'tex');
legend('x_4', 'x̂_4', 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/x4_4_lqe_plot4.png');

% График ошибки оценки для системы 4: (aQ,aR)
figure;
plot(out.xhat4 - out.x4, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('e(t)', 'Interpreter', 'tex');
legend('e_1', 'e_2', 'e_3', 'e_4', 'Location', 'southeast', 'Interpreter', 'tex');
title('Ошибка наблюдателя e(t) при (aQ,aR)', 'Interpreter', 'tex');
saveas(gcf, 'images/err4_lqe_plot4.png');
