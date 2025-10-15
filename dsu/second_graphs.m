t = out.x1.Time;
x1 = out.x1.Data;
u1 = out.u1.Data;
x2 = out.x2.Data;
u2 = out.u2.Data;
x3 = out.x3.Data;
u3 = out.u3.Data;
x4 = out.x4.Data;
u4 = out.u4.Data;
x5 = out.x5.Data;
u5 = out.u5.Data;

% График координаты x1 (Случай 1)
figure;
plot(t, x1(:,1), 'b-', 'LineWidth', 1.5);
hold on;
plot(t, x1(:,2), 'r-', 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('x(t)');
title('Случай 1: z_1 = 0.1, z_2 = 0.7');
legend('x_1', 'x_2');
saveas(gcf, 'first_case.png');

% График координаты x2 (Случай 2)
figure;
plot(t, x2(:,1), 'b-', 'LineWidth', 1.5);
hold on;
plot(t, x2(:,2), 'r-', 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('x(t)');
title('Случай 2: z_1 = -1.2, z_2 = -0.4');
legend('x_1', 'x_2');
saveas(gcf, 'second_case.png');

% График координаты x3 (Случай 3)
figure;
plot(t, x3(:,1), 'b-', 'LineWidth', 1.5);
hold on;
plot(t, x3(:,2), 'r-', 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('x(t)');
title('Случай 3: z_1 = 0.1, z_2 = 0.5');
legend('x_1', 'x_2');
saveas(gcf, 'third_case.png');

% График координаты x4 (Случай 4)
figure;
plot(t, x4(:,1), 'b-', 'LineWidth', 1.5);
hold on;
plot(t, x4(:,2), 'r-', 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('x(t)');
title('Случай 4: z_{12} = ± 1.2j');
legend('x_1', 'x_2');
saveas(gcf, 'fourth_case.png');

% График координаты x5 (Случай 5)
figure;
plot(t, x5(:,1), 'b-', 'LineWidth', 1.5);
hold on;
plot(t, x5(:,2), 'r-', 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('x(t)');
title('Случай 5: z_{12} = -0.8 ± 0.7j');
legend('x_1', 'x_2');
saveas(gcf, 'fifth_case.png');

% График управления u1 (Случай 1)
figure;
plot(t, u1, 'k-', 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('u(t)');
title('Управление - Случай 1: z_1 = 0.1, z_2 = 0.7');
saveas(gcf, 'first_case_u.png');

% График управления u2 (Случай 2)
figure;
plot(t, u2, 'k-', 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('u(t)');
title('Управление - Случай 2: z_1 = -1.2, z_2 = -0.4');
saveas(gcf, 'second_case_u.png');

% График управления u3 (Случай 3)
figure;
plot(t, u3, 'k-', 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('u(t)');
title('Управление - Случай 3: z_1 = 0.1, z_2 = 0.5');
saveas(gcf, 'third_case_u.png');

% График управления u4 (Случай 4)
figure;
plot(t, u4, 'k-', 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('u(t)');
title('Управление - Случай 4: z_{12} = ± 1.2j');
saveas(gcf, 'fourth_case_u.png');

% График управления u5 (Случай 5)
figure;
plot(t, u5, 'k-', 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('u(t)');
title('Управление - Случай 5: z_{12} = -0.8 ± 0.7j');
saveas(gcf, 'fifth_case_u.png');