time = out.x1_1.Time;
x1_1 = out.x1_1.Data;
x2_1 = out.x2_1.Data;
x1_2 = out.x1_2.Data;
x2_2 = out.x2_2.Data;
x1_3 = out.x1_3.Data;
x2_3 = out.x2_3.Data;
x3_1 = out.x3_1.Data;
x4_1 = out.x4_1.Data;
x3_2 = out.x3_2.Data;
x4_2 = out.x4_2.Data;
x3_3 = out.x3_3.Data;
x4_3 = out.x4_3.Data;

% 1

figure;
plot(time, x3_1(:,1), 'LineWidth', 2); hold on;
plot(time, x3_1(:,2), 'LineWidth', 2);
plot(time, x3_1(:,3), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
xlim([0 5]);
legend('x_1(t)','x_2(t)','x_3(t)');
title('Вектор состояния системы с регулятором K_{3,1}');
saveas(gcf, 'images/reg_x_3_1.png');

figure;
plot(time, x4_1(:,1), 'LineWidth', 2); hold on;
plot(time, x4_1(:,2), 'LineWidth', 2);
plot(time, x4_1(:,3), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_1(t)','x_2(t)','x_3(t)');
title('Вектор состояния системы с регулятором K_{4,1}');
saveas(gcf, 'images/reg_x_4_1.png');

figure;
plot(time, K3_1*x3_1', 'LineWidth', 2); hold on;
plot(time, K4_1*x4_1', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('u(t)');
legend('u(t) при K_{3,1}','u(t) при K_{4,1}');
title('Управление системой с регуляторами K_{3,1} и K_{4,1}');
saveas(gcf, 'images/reg_u_K31_K41.png');

figure;
plot(time, K1_1*x1_1', 'LineWidth', 2); hold on;
plot(time, K2_1*x2_1', 'LineWidth', 2);
plot(time, K3_1*x3_1', 'LineWidth', 2);
plot(time, K4_1*x4_1', 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t');
ylabel('u(t)');
legend('u(t) при K_{1,1}','u(t) при K_{2,1}','u(t) при K_{3,1}','u(t) при K_{4,1}');
title('Управление системой с регуляторами K_{1,1}, K_{2,1}, K_{3,1} и K_{4,1}');
saveas(gcf, 'images/reg_u_K11_K21_K31_K41.png');

% 2

figure;
plot(time, x3_2(:,1), 'LineWidth', 2); hold on;
plot(time, x3_2(:,2), 'LineWidth', 2);
plot(time, x3_2(:,3), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
xlim([0 3]);
legend('x_1(t)','x_2(t)','x_3(t)');
title('Вектор состояния системы с регулятором K_{3,2}');
saveas(gcf, 'images/reg_x_3_2.png');

figure;
plot(time, x4_2(:,1), 'LineWidth', 2); hold on;
plot(time, x4_2(:,2), 'LineWidth', 2);
plot(time, x4_2(:,3), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_1(t)','x_2(t)','x_3(t)');
xlim([0 5]);
title('Вектор состояния системы с регулятором K_{4,2}');
saveas(gcf, 'images/reg_x_4_2.png');

figure;
plot(time, K3_2*x3_2', 'LineWidth', 2); hold on;
plot(time, K4_2*x4_2', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('u(t)');
xlim([0 5]);
legend('u(t) при K_{3,2}','u(t) при K_{4,2}');
title('Управление системой с регуляторами K_{3,2} и K_{4,2}');
saveas(gcf, 'images/reg_u_K32_K42.png');

figure;
plot(time, K1_2*x1_2', 'LineWidth', 2); hold on;
plot(time, K2_2*x2_2', 'LineWidth', 2);
plot(time, K3_2*x3_2', 'LineWidth', 2);
plot(time, K4_2*x4_2', 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t');
ylabel('u(t)');
xlim([0 5]);
legend('u(t) при K_{1,2}','u(t) при K_{2,2}','u(t) при K_{3,2}','u(t) при K_{4,2}');
title('Управление системой с регуляторами K_{1,2}, K_{2,2}, K_{3,2} и K_{4,2}');
saveas(gcf, 'images/reg_u_K12_K22_K32_K42.png');

% 3

figure;
plot(time, x3_3(:,1), 'LineWidth', 2); hold on;
plot(time, x3_3(:,2), 'LineWidth', 2);
plot(time, x3_3(:,3), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
xlim([0 3]);
legend('x_1(t)','x_2(t)','x_3(t)');
title('Вектор состояния системы с регулятором K_{3,3}');
saveas(gcf, 'images/reg_x_3_3.png');

figure;
plot(time, x4_3(:,1), 'LineWidth', 2); hold on;
plot(time, x4_3(:,2), 'LineWidth', 2);
plot(time, x4_3(:,3), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_1(t)','x_2(t)','x_3(t)');
xlim([0 3]);
title('Вектор состояния системы с регулятором K_{4,3}');
saveas(gcf, 'images/reg_x_4_3.png');

figure;
plot(time, K3_3*x3_3', 'LineWidth', 2); hold on;
plot(time, K4_3*x4_3', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('u(t)');
xlim([0 3]);
legend('u(t) при K_{3,3}','u(t) при K_{4,3}');
title('Управление системой с регуляторами K_{3,3} и K_{4,3}');
saveas(gcf, 'images/reg_u_K33_K43.png');

figure;
plot(time, K1_3*x1_3', 'LineWidth', 2); hold on;
plot(time, K2_3*x2_3', 'LineWidth', 2);
plot(time, K3_3*x3_3', 'LineWidth', 2);
plot(time, K4_3*x4_3', 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t');
ylabel('u(t)');
xlim([0 3]);
legend('u(t) при K_{1,3}','u(t) при K_{2,3}','u(t) при K_{3,3}','u(t) при K_{4,3}');
title('Управление системой с регуляторами K_{1,3}, K_{2,3}, K_{3,3} и K_{4,3}');
saveas(gcf, 'images/reg_u_K13_K23_K33_K43.png');