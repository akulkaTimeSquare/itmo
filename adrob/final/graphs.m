%% full
close all;
figure;
plot(out.u, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('u(t)', 'Interpreter', 'tex');
title('Управление u(t) при градиентном АА', 'Interpreter', 'tex');
legend('u', 'Interpreter', 'tex');
saveas(gcf, 'images/u_aa.png');

figure;
plot(out.y, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('y(t)', 'Interpreter', 'tex');
title('Выходная переменная y(t) при градиентном АА', 'Interpreter', 'tex');
legend('e_f', 'Interpreter', 'tex');
saveas(gcf, 'images/y_aa.png');

figure;
plot(out.ef, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('e_f(t)', 'Interpreter', 'tex');
title('Ошибка компенсации воздействия f(t) при градиентном АА', 'Interpreter', 'tex');
legend('e_f', 'Interpreter', 'tex');
saveas(gcf, 'images/ef_aa.png');

figure;
plot(out.x, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояния x(t) при градиентном АА', 'Interpreter', 'tex');
legend('x_1', 'x_2', 'Interpreter', 'tex');
saveas(gcf, 'images/x_aa.png');

figure;
plot(out.ksi, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('\xi(t)', 'Interpreter', 'tex');
title('Наблюдатель состояния \xi_f(t) при градиентном АА', 'Interpreter', 'tex');
legend('\xi_1', '\xi_2', 'Interpreter', 'tex');
saveas(gcf, 'images/xi_aa.png');

figure;
hold on;
plot(out.theta_hat, 'LineWidth', 2);
title('Оценка параметров \theta при градиентном АА');
xlabel('t'); ylabel('\theta');
legend('\theta_1', '\theta_2', '\theta_3', '\theta_4', 'Interpreter', 'tex', 'Location', 'southwest');
grid on;
saveas(gcf, 'images/theta_aa.png');

%% t = 25
close all;
figure;
plot(out.u, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('u(t)', 'Interpreter', 'tex');
title('Управление u(t) при градиентном АА', 'Interpreter', 'tex');
legend('u', 'Interpreter', 'tex');
saveas(gcf, 'images/u_aa_25.png');

figure;
plot(out.y, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('y(t)', 'Interpreter', 'tex');
title('Выходная переменная y(t) при градиентном АА', 'Interpreter', 'tex');
legend('e_f', 'Interpreter', 'tex');
saveas(gcf, 'images/y_aa_25.png');

figure;
plot(out.ef, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('e_f(t)', 'Interpreter', 'tex');
title('Ошибка компенсации воздействия f(t) при градиентном АА', 'Interpreter', 'tex');
legend('e_f', 'Interpreter', 'tex');
saveas(gcf, 'images/ef_aa_25.png');

figure;
plot(out.x, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояния x(t) при градиентном АА', 'Interpreter', 'tex');
legend('x_1', 'x_2', 'Interpreter', 'tex');
saveas(gcf, 'images/x_aa_25.png');

figure;
plot(out.ksi, 'LineWidth', 2);
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('\xi(t)', 'Interpreter', 'tex');
title('Наблюдатель состояния \xi_f(t)', 'Interpreter', 'tex');
legend('\xi_1', '\xi_2', 'Interpreter', 'tex');
saveas(gcf, 'images/xi_aa_25.png');

figure;
hold on;
plot(out.theta_hat, 'LineWidth', 2);
title('Оценка параметров \theta при градиентном АА');
xlabel('t'); ylabel('\theta');
grid on;
legend('\theta_1', '\theta_2', '\theta_3', '\theta_4', 'Interpreter', 'tex');
saveas(gcf, 'images/theta_aa_25.png');