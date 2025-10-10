time = out.x.Time;
x1 = out.x1.Data;
x2 = out.x2.Data;
x3 = out.x3.Data;
hatx1 = out.hatx1.Data;
hatx2 = out.hatx2.Data;
hatx3 = out.hatx3.Data;
u1 = out.u1.Data;
u2 = out.u2.Data;
u3 = out.u3.Data;

%1

figure;
plot(time, x1(:,1), 'LineWidth', 2); hold on;
plot(time, hatx1(:,1), 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t');
ylabel('x_1(t)');
legend('x_1(t)','x̂_1(t)');
xlim([0 1]);
title('Вторая компонента векторов состояний при \alpha_{K1} = \alpha_{L1} = 12', 'Interpreter', 'tex');
saveas(gcf, 'images/obsreg_x1_1.png');

figure;
plot(time, x1(:,2), 'LineWidth', 2); hold on;
plot(time, hatx1(:,2), 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t');
ylabel('x_2(t)');
legend('x_2(t)','x̂_2(t)');
xlim([0 1]);
title('Вторая компонента векторов состояний при \alpha_{K1} = \alpha_{L1} = 12', 'Interpreter', 'tex');
saveas(gcf, 'images/obsreg_x2_1.png');

figure;
plot(time, x1(:,3), 'LineWidth', 2); hold on;
plot(time, hatx1(:,3), 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t');
ylabel('x_3(t)');
legend('x_3(t)','x̂_3(t)');
xlim([0 1]);
title('Третья компонента векторов состояний при \alpha_{K1} = \alpha_{L1} = 12', 'Interpreter', 'tex');
saveas(gcf, 'images/obsreg_x3_1.png');

figure;
plot(time, x1(:,4), 'LineWidth', 2); hold on;
plot(time, hatx1(:,4), 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t');
ylabel('x_4(t)');
legend('x_4(t)','x̂_4(t)');
xlim([0 1]);
title('Четвертая компонента векторов состояний при \alpha_{K1} = \alpha_{L1} = 12', 'Interpreter', 'tex');
saveas(gcf, 'images/obsreg_x4_1.png');

figure;
plot(time, u1, 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('u(t)');
xlim([0 1]);
title('Управление системой при \alpha_{K1} = \alpha_{L1} = 12', 'Interpreter', 'tex');
saveas(gcf, 'images/obsreg_u_1.png');

figure;
plot(time, x1(:,1) - hatx1(:,1), 'LineWidth', 2); hold on;
plot(time, x1(:,2) - hatx1(:,2), 'LineWidth', 2);
plot(time, x1(:,3) - hatx1(:,3), 'LineWidth', 2);
plot(time, x1(:,4) - hatx1(:,4), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('e(t)');
legend('e_1(t)','e_2(t)','e_3(t)','e_4(t)');
xlim([0 1]);
title('Ошибка наблюдателя e(t) = x(t) - x̂(t) при \alpha_{K1} = \alpha_{L1} = 12', 'Interpreter', 'tex');
saveas(gcf, 'images/obsreg_e_1.png');

%2

figure;
plot(time, x2(:,1), 'LineWidth', 2); hold on;
plot(time, hatx2(:,1), 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t');
ylabel('x_1(t)');
legend('x_1(t)','x̂_1(t)');
xlim([0 5]);
title('Вторая компонента векторов состояний при \alpha_{K2} = 12 и \alpha_{L2} = 1', 'Interpreter', 'tex');
saveas(gcf, 'images/obsreg_x1_2.png');

figure;
plot(time, x2(:,2), 'LineWidth', 2); hold on;
plot(time, hatx2(:,2), 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t');
ylabel('x_2(t)');
legend('x_2(t)','x̂_2(t)');
xlim([0 5]);
title('Вторая компонента векторов состояний при \alpha_{K2} = 12 и \alpha_{L2} = 1', 'Interpreter', 'tex');
saveas(gcf, 'images/obsreg_x2_2.png');

figure;
plot(time, x2(:,3), 'LineWidth', 2); hold on;
plot(time, hatx2(:,3), 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t');
ylabel('x_3(t)');
legend('x_3(t)','x̂_3(t)');
xlim([0 5]);
title('Третья компонента векторов состояний при \alpha_{K2} = 12 и \alpha_{L2} = 1', 'Interpreter', 'tex');
saveas(gcf, 'images/obsreg_x3_2.png');

figure;
plot(time, x2(:,4), 'LineWidth', 2); hold on;
plot(time, hatx2(:,4), 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t');
ylabel('x_4(t)');
legend('x_4(t)','x̂_4(t)');
xlim([0 5]);
title('Четвертая компонента векторов состояний при \alpha_{K2} = 12 и \alpha_{L2} = 1', 'Interpreter', 'tex');
saveas(gcf, 'images/obsreg_x4_2.png');

figure;
plot(time, u2, 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('u(t)');
xlim([0 5]);
title('Управление системой при \alpha_{K2} = 12 и \alpha_{L2} = 1', 'Interpreter', 'tex');
saveas(gcf, 'images/obsreg_u_2.png');

figure;
plot(time, x2(:,1) - hatx2(:,1), 'LineWidth', 2); hold on;
plot(time, x2(:,2) - hatx2(:,2), 'LineWidth', 2);
plot(time, x2(:,3) - hatx2(:,3), 'LineWidth', 2);
plot(time, x2(:,4) - hatx2(:,4), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('e(t)');
xlim([0 5]);
legend('e_1(t)','e_2(t)','e_3(t)','e_4(t)');
title('Ошибка наблюдателя e(t) = x(t) - x̂(t) при \alpha_{K2} = 12 и \alpha_{L2} = 1', 'Interpreter', 'tex');
saveas(gcf, 'images/obsreg_e_2.png');

%3

figure;
plot(time, x3(:,1), 'LineWidth', 2); hold on;
plot(time, hatx3(:,1), 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t');
ylabel('x_1(t)');
legend('x_1(t)','x̂_1(t)');
xlim([0 5]);
title('Вторая компонента векторов состояний при \alpha_{K3} = 1 и \alpha_{L3} = 12', 'Interpreter', 'tex');
saveas(gcf, 'images/obsreg_x1_3.png');

figure;
plot(time, x3(:,2), 'LineWidth', 2); hold on;
plot(time, hatx3(:,2), 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t');
ylabel('x_2(t)');
legend('x_2(t)','x̂_2(t)');
xlim([0 5]);
title('Вторая компонента векторов состояний при \alpha_{K3} = 1 и \alpha_{L3} = 12', 'Interpreter', 'tex');
saveas(gcf, 'images/obsreg_x2_3.png');

figure;
plot(time, x3(:,3), 'LineWidth', 2); hold on;
plot(time, hatx3(:,3), 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t');
ylabel('x_3(t)');
legend('x_3(t)','x̂_3(t)');
xlim([0 5]);
title('Третья компонента векторов состояний при \alpha_{K3} = 1 и \alpha_{L3} = 12', 'Interpreter', 'tex');
saveas(gcf, 'images/obsreg_x3_3.png');  

figure;
plot(time, x3(:,4), 'LineWidth', 2); hold on;
plot(time, hatx3(:,4), 'LineWidth', 2, 'LineStyle', '--');
grid on;
xlabel('t');
ylabel('x_4(t)');
legend('x_4(t)','x̂_4(t)');
xlim([0 5]);
title('Четвертая компонента векторов состояний при \alpha_{K3} = 1 и \alpha_{L3} = 12', 'Interpreter', 'tex');
saveas(gcf, 'images/obsreg_x4_3.png');

figure;
plot(time, u3, 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('u(t)');
xlim([0 5]);
title('Управление системой при \alpha_{K3} = 1 и \alpha_{L3} = 12', 'Interpreter', 'tex');
saveas(gcf, 'images/obsreg_u_3.png');

figure;
plot(time, x3(:,1) - hatx3(:,1), 'LineWidth', 2); hold on;
plot(time, x3(:,2) - hatx3(:,2), 'LineWidth', 2);
plot(time, x3(:,3) - hatx3(:,3), 'LineWidth', 2);
plot(time, x3(:,4) - hatx3(:,4), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('e(t)');
legend('e_1(t)','e_2(t)','e_3(t)','e_4(t)');
xlim([0 1]);
title('Ошибка наблюдателя e(t) = x(t) - x̂(t) при \alpha_{K3} = 1 и \alpha_{L3} = 12', 'Interpreter', 'tex');
saveas(gcf, 'images/obsreg_e_3.png');

figure;
plot(time, u1, 'LineWidth', 2); hold on;
plot(time, u2, 'LineWidth', 2);
plot(time, u3, 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('u(t)');
legend('u(t) при \alpha_K = \alpha_L', 'u(t) при \alpha_K > \alpha_L','u(t) при \alpha_K < \alpha_L');
xlim([0 3]);
title('Управление системой при различных спектрах');
saveas(gcf, 'images/obsreg_u_all.png');