time = out.y_hat.Time;
y_hat = out.y_hat.Data;
e = out.e.Data;
theta_tilde = out.theta_tilde.Data;
e_x = out.e_x.Data;
norm_e_x = sqrt(e_x(:, 1).^2 + e_x(:, 2).^2);

figure;
plot(time, norm_e_x, 'LineWidth', 2); 
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('||e_x(t)||', 'Interpreter', 'tex');
title('Норма ошибки оценки переменных состояния', 'Interpreter', 'tex');
saveas(gcf, 'images/gamma3_ex.png');

figure;
plot(time, theta_tilde(:, 1), 'LineWidth', 2); hold on;
plot(time, theta_tilde(:, 2), 'LineWidth', 2);
plot(time, theta_tilde(:, 3), 'LineWidth', 2);
plot(time, theta_tilde(:, 4), 'LineWidth', 2);
hold off;
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('\theta(t)', 'Interpreter', 'tex');
title('Параметрические ошибки оценки \theta', 'Interpreter', 'tex');
legend("\theta_1", "\theta_2", "\theta_3", "\theta_4", 'Location', 'northeast', 'Interpreter', 'tex');
saveas(gcf, 'images/gamma3_theta.png');