t = out.x.Time;
x0 = out.x1.Data;
x1 = out.x1.Data;
x2 = out.x2.Data;
x3 = out.x3.Data;

%% График x1(t) для всех систем на одном рисунке
figure;
plot(t, x0(:, 1), 'LineWidth', 2, "DisplayName", "(r_0, k_0)"); hold on;
plot(t, x1(:, 1), 'LineWidth', 2, 'LineStyle', "--", "DisplayName", "(r_1, k_1)");
plot(t, x2(:, 1), 'LineWidth', 2, "DisplayName", "(r_2, k_2)");
plot(t, x3(:, 1), 'LineWidth', 2, "DisplayName", "(r_3, k_3)");
hold off;
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояния x_1 замкнутых систем', 'Interpreter', 'tex');
legend('Location', 'best', 'Interpreter', 'tex');
saveas(gcf, 'images/x1_m.png');

%% График x2(t) для всех систем на одном рисунке
figure;
plot(t, x0(:, 2), 'LineWidth', 2, "DisplayName", "(r_0, k_0)"); hold on;
plot(t, x1(:, 2), 'LineWidth', 2, 'LineStyle', "--", "DisplayName", "(r_1, k_1)");
plot(t, x2(:, 2), 'LineWidth', 2, "DisplayName", "(r_2, k_2)");
plot(t, x3(:, 2), 'LineWidth', 2, "DisplayName", "(r_3, k_3)");
hold off;
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('x(t)', 'Interpreter', 'tex');
title('Состояния x_2 объекта замкнутых систем', 'Interpreter', 'tex');
legend('Location', 'best', 'Interpreter', 'tex');
saveas(gcf, 'images/x2_m.png');

%% График всех управлений u(t) на одном рисунке
figure;
plot(out.u, 'LineWidth', 2, 'DisplayName', '(r_0, k_0)'); hold on;
plot(out.u1, 'LineWidth', 2, 'DisplayName', '(r_1, k_1)', 'LineStyle', '--');
plot(out.u2, 'LineWidth', 2, 'DisplayName', '(r_2, k_2)');
plot(out.u3, 'LineWidth', 2, 'DisplayName', '(r_3, k_3)');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('u(t)', 'Interpreter', 'tex');
title('Оптимальные управления замкнутых систем', 'Interpreter', 'tex');
legend('Location', 'best', 'Interpreter', 'tex');
saveas(gcf, 'images/u_m.png');