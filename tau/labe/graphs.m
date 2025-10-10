time = out.xw.Time;
xw = out.xw.Data;
xwhat = out.xwhat.Data;
u = out.u.Data';

A = [0 1; -1 1];
B = [1 2; 1 0];
C = [1 -2; 0 3];
D = [-3 0; 0 1];

Bf = [1 2; 1 3];
Df = [1 0; 0 1];
Cz = [1 2; 4 0];
Dz = [4 0; 0 1];

Yg = [0 0 2 0 0 0 0 0;
      0 0 0 -3 0 0 0 0];
Y1 = [0 0 0 0 0 0 7 0;
      0 -3 0 0 0 0 0 0];
Y2 = [0 0 0 0 0 -5 0 0;
      0 0 0 0 0 0 3 0];

w = xw(:, 1:8)';
what = xwhat(:, 1:8)';

x = xw(:, 9:10)';
xhat = xwhat(:, 9:10)';

g = Yg*w;
f1 = Y1*w;
f2 = Y2*w;
z = Cz*x + Dz*u - g;
y = C*x + D*u + Df*f2;

% График 1: Состояния системы и их оценки
figure;
% Состояния системы (сплошные линии)
plot(time, x(1, :), 'LineWidth', 2);
hold on;
plot(time, x(2, :), 'LineWidth', 2);
plot(time, xhat(1, :), 'LineWidth', 2, 'LineStyle', '--');
plot(time, xhat(2, :), 'LineWidth', 2, 'LineStyle', '--');
% Оценки состояния (пунктирные линии)
legend('x_1', 'x_2', 'x̂_1', 'x̂_2', 'Location', 'best');
title('Состояния системы и их оценки');
xlabel('t');
ylabel('x(t)');
grid on;
% Добавляем отступы по Y
yl = ylim;
ylim([yl(1) - 0.05*diff(yl), yl(2) + 0.05*diff(yl)]);
saveas(gcf, 'images/x.png');

figure;
plot(time, w(1, :), 'LineWidth', 2);
hold on;
plot(time, w(2, :), 'LineWidth', 2);
plot(time, what(1, :), 'LineWidth', 2, 'LineStyle', '--');
plot(time, what(2, :), 'LineWidth', 2, 'LineStyle', '--');
legend('w_1', 'w_2', 'ŵ_1', 'ŵ_2', 'Location', 'southeast');
title('Первая и вторая компоненты w(t) и их оценки ŵ(t)');
xlabel('t');
ylabel('w(t)');
grid on;
% Добавляем отступы по Y
yl = ylim;
ylim([yl(1) - 0.05*diff(yl), yl(2) + 0.05*diff(yl)]);
saveas(gcf, 'images/w12.png');

figure;
plot(time, w(3, :), 'LineWidth', 2);
hold on;
plot(time, w(4, :), 'LineWidth', 2);
plot(time, what(3, :), 'LineWidth', 2, 'LineStyle', '--');
plot(time, what(4, :), 'LineWidth', 2, 'LineStyle', '--');
legend('w_3', 'w_4', 'ŵ_3', 'ŵ_4', 'Location', 'best');
title('Третья и четвертая компоненты w(t) и их оценки ŵ(t)');
xlabel('t');
ylabel('w(t)');
grid on;
% Добавляем отступы по Y
yl = ylim;
ylim([yl(1) - 0.05*diff(yl), yl(2) + 0.05*diff(yl)]);
saveas(gcf, 'images/w34.png');

figure;
plot(time, w(5, :), 'LineWidth', 2);
hold on;
plot(time, w(6, :), 'LineWidth', 2);
plot(time, what(5, :), 'LineWidth', 2, 'LineStyle', '--');
plot(time, what(6, :), 'LineWidth', 2, 'LineStyle', '--');
legend('w_5', 'w_6', 'ŵ_5', 'ŵ_6', 'Location', 'southeast');
title('Пятая и шестая компоненты w(t) и их оценки ŵ(t)');
xlabel('t');
ylabel('w(t)');
grid on;
% Добавляем отступы по Y
yl = ylim;
ylim([yl(1) - 0.05*diff(yl), yl(2) + 0.05*diff(yl)]);
saveas(gcf, 'images/w56.png');

figure;
plot(time, w(7, :), 'LineWidth', 2);
hold on;
plot(time, w(8, :), 'LineWidth', 2);
plot(time, what(7, :), 'LineWidth', 2, 'LineStyle', '--');
plot(time, what(8, :), 'LineWidth', 2, 'LineStyle', '--');
legend('w_7', 'w_8', 'ŵ_7', 'ŵ_8', 'Location', 'northeast');
title('Седьмая и восьмая компоненты w(t) и их оценки ŵ(t)');
xlabel('t');
ylabel('w(t)');
grid on;
% Добавляем отступы по Y
yl = ylim;
ylim([yl(1) - 0.05*diff(yl), yl(2) + 0.05*diff(yl)]);
saveas(gcf, 'images/w78.png');

figure;
plot(time, f1(1, :), 'LineWidth', 2);
hold on;
plot(time, f1(2, :), 'LineWidth', 2);
legend('f_1_1', 'f_1_2', 'Location', 'best');
title('Воздействие f_1(t) на состояния системы');
xlabel('t');
ylabel('f_1(t)');
grid on;
% Добавляем отступы по Y
yl = ylim;
ylim([yl(1) - 0.05*diff(yl), yl(2) + 0.05*diff(yl)]);
saveas(gcf, 'images/f1.png');

figure;
plot(time, f2(1, :), 'LineWidth', 2);
hold on;
plot(time, f2(2, :), 'LineWidth', 2);
legend('f_2_1', 'f_2_2', 'Location', 'northeast');
title('Воздействие f_2(t) на выход системы');
xlabel('t');
ylabel('f_2(t)');
grid on;
% Добавляем отступы по Y
yl = ylim;
ylim([yl(1) - 0.05*diff(yl), yl(2) + 0.05*diff(yl)]);
saveas(gcf, 'images/f2.png');

figure;
plot(time, g(1, :), 'LineWidth', 2);
hold on;
plot(time, g(2, :), 'LineWidth', 2);
legend('g_1', 'g_2', 'Location', 'best');
title('Задающий сигнал g(t) системы');
xlabel('t');
ylabel('g(t)');
grid on;
% Добавляем отступы по Y
yl = ylim;
ylim([yl(1) - 0.05*diff(yl), yl(2) + 0.05*diff(yl)]);
saveas(gcf, 'images/g.png');

% График 2: Управление u
figure;
plot(time, u(1, :), 'LineWidth', 2);
hold on;
plot(time, u(2, :), 'LineWidth', 2);
legend('u_1', 'u_2', 'Location', 'best');
title('Управление многоканальной системы u(t)');
xlabel('t');
ylabel('u(t)');
grid on;
% Добавляем отступы по Y
yl = ylim;
ylim([yl(1) - 0.05*diff(yl), yl(2) + 0.05*diff(yl)]);
saveas(gcf, 'images/u.png');

% График 3: Виртуальный выход z
figure;
plot(time, z(1, :), 'LineWidth', 2);
hold on;
plot(time, z(2, :), 'LineWidth', 2);
legend('z_1', 'z_2', 'Location', 'best');
title('Виртуальный выход многоканальной системы z(t)');
xlabel('t');
ylabel('z(t)');
grid on;
% Добавляем отступы по Y
yl = ylim;
ylim([yl(1) - 0.05*diff(yl), yl(2) + 0.05*diff(yl)]);
saveas(gcf, 'images/z.png');

figure;
plot(time, y(1, :), 'LineWidth', 2);
hold on;
plot(time, y(2, :), 'LineWidth', 2);
legend('y_1', 'y_2', 'Location', 'best');
title('Фактический выход многоканальной системы y(t)');
xlabel('t');
ylabel('y(t)');
grid on;
% Добавляем отступы по Y
yl = ylim;
ylim([yl(1) - 0.05*diff(yl), yl(2) + 0.05*diff(yl)]);
saveas(gcf, 'images/y.png');

figure;
plot(time, x(1, :) - xhat(1, :), 'LineWidth', 2);
hold on;
plot(time, x(2, :) - xhat(2, :), 'LineWidth', 2);
legend('e_x_1', 'e_x_2', 'Location', 'northeast');
title('Ошибка оценивания состояния системы e_x(t)');
xlabel('t');
ylabel('e_x(t)');
grid on;
% Добавляем отступы по Y
yl = ylim;
ylim([yl(1) - 0.05*diff(yl), yl(2) + 0.05*diff(yl)]);
saveas(gcf, 'images/ex.png');

figure;
plot(time, w(1, :) - what(1, :), 'LineWidth', 2);
hold on;
plot(time, w(2, :) - what(2, :), 'LineWidth', 2);
plot(time, w(3, :) - what(3, :), 'LineWidth', 2);
plot(time, w(4, :) - what(4, :), 'LineWidth', 2);
legend('e_w_1', 'e_w_2', 'e_w_3', 'e_w_4', 'Location', 'northeast');
title('Ошибка оценивания первых четырех воздействий e_w(t)');
xlabel('t');
ylabel('e_w_i(t)');
grid on;
% Добавляем отступы по Y
yl = ylim;
ylim([yl(1) - 0.05*diff(yl), yl(2) + 0.05*diff(yl)]);
saveas(gcf, 'images/ew14.png');

figure;
plot(time, w(5, :) - what(5, :), 'LineWidth', 2);
hold on;
plot(time, w(6, :) - what(6, :), 'LineWidth', 2);
plot(time, w(7, :) - what(7, :), 'LineWidth', 2);
plot(time, w(8, :) - what(8, :), 'LineWidth', 2);
legend('e_w_5', 'e_w_6', 'e_w_7', 'e_w_8', 'Location', 'northeast');
title('Ошибка оценивания последних четырех воздействий e_w(t)');
xlabel('t');
ylabel('e_w_i(t)');
grid on;
% Добавляем отступы по Y
yl = ylim;
ylim([yl(1) - 0.05*diff(yl), yl(2) + 0.05*diff(yl)]);
saveas(gcf, 'images/ew58.png');
