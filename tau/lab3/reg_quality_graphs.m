time = out.x1.Time;
x1 = out.x1.Data;
x2 = out.x2.Data;
x3 = out.x3.Data;
x4 = out.x4.Data;

figure;
plot(time, x1(:,1), 'LineWidth', 2); hold on;
plot(time, x1(:,2), 'LineWidth', 2);
plot(time, x1(:,3), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_1(t)','x_2(t)','x_3(t)');
title('Вектор состояния системы при Q = I, R = 1');
saveas(gcf, 'images/reg_quality_x_1.png');

figure;
plot(time, K1*x1', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('u(t)');
title('Управление системой при Q = I, R = 1');
saveas(gcf, 'images/reg_quality_u_1.png');

figure;
plot(time, x2(:,1), 'LineWidth', 2); hold on;
plot(time, x2(:,2), 'LineWidth', 2);
plot(time, x2(:,3), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_1(t)','x_2(t)','x_3(t)');
title('Вектор состояния системы при Q = I, R = 0');
saveas(gcf, 'images/reg_quality_x_2.png');

figure;
plot(time, K2*x2', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('u(t)');
title('Управление системой при Q = I, R = 0');
saveas(gcf, 'images/reg_quality_u_2.png');

figure;
plot(time, x3(:,1), 'LineWidth', 2); hold on;
plot(time, x3(:,2), 'LineWidth', 2);
plot(time, x3(:,3), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_1(t)','x_2(t)','x_3(t)');
title('Вектор состояния системы при Q = 0, R = 1');
saveas(gcf, 'images/reg_quality_x_3.png');

figure;
plot(time, K3*x3', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('u(t)');
title('Управление системой при Q = 0, R = 1');
saveas(gcf, 'images/reg_quality_u_3.png');

figure;
plot(time, x4(:,1), 'LineWidth', 2); hold on;
plot(time, x4(:,2), 'LineWidth', 2);
plot(time, x4(:,3), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_1(t)','x_2(t)','x_3(t)');
title('Вектор состояния системы при Q = 0, R = 0');
saveas(gcf, 'images/reg_quality_x_4.png');

figure;
plot(time, K4*x4', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('u(t)');
title('Управление системой при Q = 0, R = 0');
saveas(gcf, 'images/reg_quality_u_4.png');

figure;
plot(time, K1*x1', 'LineWidth', 2); hold on;
plot(time, K2*x2', 'LineWidth', 2);
plot(time, K3*x3', 'LineWidth', 2);
plot(time, K4*x4', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('u(t)');
legend('u(t) при (Q = I, R = 1)','u(t) при (Q = I, R = 0)','u(t) при (Q = 0, R = 1)','u(t) при (Q = 0, R = 0)', 'Location', 'southeast');
title('Управление регуляторами при различных парах (Q, R)');
saveas(gcf, 'images/reg_quality_u_K1_K2_K3_K4.png');

beta = -3;
r = 2;
figure;
th = 0:pi/50:2*pi;
xunit = r * cos(th) + beta;
yunit = r * sin(th);
plot(xunit, yunit, "LineWidth", 2);
hold on
plot(real(e1), imag(e1), "o", "LineWidth", 8); hold on;
plot(real(e2), imag(e2), "x", "LineWidth", 6);
plot(real(e3), imag(e3), "+", "LineWidth", 4);
plot(real(e4), imag(e4), "s", "LineWidth", 2);
axis equal
grid on
xlabel("Re")
ylabel("Im")
legend("Круг r, (\beta, 0)", "Q = I, R = 1", "Q = I, R = 0", "Q = 0, R = 1", "Q = 0, R = 0", "Location", "southeast")
title("Комплексная плоскость и собственные числа (все случаи)")
hold off
saveas(gcf, "images/reg_quality_all.png")