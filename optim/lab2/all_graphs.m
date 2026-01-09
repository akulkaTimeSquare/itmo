%% Графики x1
figure;
plot(t_vals, x1_num, 'LineWidth', 2); hold on;
plot(t_vals1, x1_num1, 'LineWidth', 2);
plot(t_vals2, x1_num2, 'LineWidth', 2);
hold off;
title('Переменные состояния x_1');
xlabel('t');
ylabel('x(t)');
grid on;
legend("x_1(t) при u^*", "x_1(t) при u_1", "x_1(t) при u_2");
saveas(gcf, "images/x1m.png");

%% Графики x2
figure;
plot(t_vals, x2_num, 'LineWidth', 2); hold on;
plot(t_vals1, x2_num1, 'LineWidth', 2);
plot(t_vals2, x2_num2, 'LineWidth', 2);
hold off;
title('Переменные состояния x_2');
xlabel('t');
ylabel('x(t)');
grid on;
legend("x_2(t) при u^*", "x_2(t) при u_1", "x_2(t) при u_2");
saveas(gcf, "images/x2m.png");

%% Графики u
figure;
plot(t_vals, u_num, 'LineWidth', 2); hold on;
plot(t_vals1, u_num1, 'LineWidth', 2);
plot(t_vals2, u_num2, 'LineWidth', 2);
hold off;
title('Оптимальные управления');
xlabel('t');
ylabel('u(t)');
legend("u^*(t)", "x_2(t) при u_1(t)", "x_2(t) при u_2(t)");
grid on;
saveas(gcf, "images/um.png");

%% Графики J
figure;
plot(t_vals, J_t_vals, 'LineWidth', 2); hold on;
plot(t_vals1, J_t_vals1, 'LineWidth', 2);
plot(t_vals2, J_t_vals2, 'LineWidth', 2);
hold off;
title('Накопленные критерии качества', 'Interpreter', 'latex');
xlabel('t');
ylabel('J(t)');
grid on;
yticks([0 1000 2000 3000 4000 5000 6000 6829])
saveas(gcf, "images/Jm.png");
