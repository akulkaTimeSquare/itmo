time = out.x1_1.Time;
x1_1 = out.x1_1.Data;
x2_1 = out.x2_1.Data;
x1_2 = out.x1_2.Data;
x2_2 = out.x2_2.Data;
x1_3 = out.x1_3.Data;
x2_3 = out.x2_3.Data;

% 1

figure;
plot(time, x1_1(:,1), 'LineWidth', 2); hold on;
plot(time, x1_1(:,2), 'LineWidth', 2);
plot(time, x1_1(:,3), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
xlim([0 3]);
legend('x_1(t)','x_2(t)','x_3(t)');
title('Вектор состояния системы с регулятором K_{1,1}');
saveas(gcf, 'images/reg_x_1_1.png');

figure;
plot(time, x2_1(:,1), 'LineWidth', 2); hold on;
plot(time, x2_1(:,2), 'LineWidth', 2);
plot(time, x2_1(:,3), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_1(t)','x_2(t)','x_3(t)');
title('Вектор состояния системы с регулятором K_{2,1}');
saveas(gcf, 'images/reg_x_2_1.png');

figure;
plot(time, K1_1*x1_1', 'LineWidth', 2); hold on;
plot(time, K2_1*x2_1', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('u(t)');
legend('u(t) при K_{1,1}','u(t) при K_{2,1}');
title('Управление системой с регуляторами K_{1,1} и K_{2,1}');
saveas(gcf, 'images/reg_u_K11_K21.png');

% 2

figure;
plot(time, x1_2(:,1), 'LineWidth', 2); hold on;
plot(time, x1_2(:,2), 'LineWidth', 2);
plot(time, x1_2(:,3), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
xlim([0 3]);
legend('x_1(t)','x_2(t)','x_3(t)');
title('Вектор состояния системы с регулятором K_{1,2}');
saveas(gcf, 'images/reg_x_1_2.png');

figure;
plot(time, x2_2(:,1), 'LineWidth', 2); hold on;
plot(time, x2_2(:,2), 'LineWidth', 2);
plot(time, x2_2(:,3), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_1(t)','x_2(t)','x_3(t)');
xlim([0 5]);
title('Вектор состояния системы с регулятором K_{2,2}');
saveas(gcf, 'images/reg_x_2_2.png');

figure;
plot(time, K1_2*x1_2', 'LineWidth', 2); hold on;
plot(time, K2_2*x2_2', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('u(t)');
xlim([0 5]);
legend('u(t) при K_{1,2}','u(t) при K_{2,2}');
title('Управление системой с регуляторами K_{1,2} и K_{2,2}');
saveas(gcf, 'images/reg_u_K12_K22.png');

% 3

figure;
plot(time, x1_3(:,1), 'LineWidth', 2); hold on;
plot(time, x1_3(:,2), 'LineWidth', 2);
plot(time, x1_3(:,3), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
xlim([0 3]);
legend('x_1(t)','x_2(t)','x_3(t)');
title('Вектор состояния системы с регулятором K_{1,3}');
saveas(gcf, 'images/reg_x_1_3.png');

figure;
plot(time, x2_3(:,1), 'LineWidth', 2); hold on;
plot(time, x2_3(:,2), 'LineWidth', 2);
plot(time, x2_3(:,3), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_1(t)','x_2(t)','x_3(t)');
xlim([0 3]);
title('Вектор состояния системы с регулятором K_{2,3}');
saveas(gcf, 'images/reg_x_2_3.png');

figure;
plot(time, K1_3*x1_3', 'LineWidth', 2); hold on;
plot(time, K2_3*x2_3', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('u(t)');
xlim([0 3]);
legend('u(t) при K_{1,3}','u(t) при K_{2,3}');
title('Управление системой с регуляторами K_{1,3} и K_{2,3}');
saveas(gcf, 'images/reg_u_K13_K23.png');