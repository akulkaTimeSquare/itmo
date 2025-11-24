load('ident_lab2_v05.mat');
% исходные параметры
a1 = zad2.a1;
a2 = zad2.a2;
b = zad2.b;
w = zad2.w;

% Параметры
Td = 0.1;
t_end = 12;
t = 0:Td:t_end;
N = t_end/Td+1;

% Входной сигнал
u1 = sin(w*t);
u2 = sin(w*t) + 0.2*sin(0.5*w*t);

% Создание дискретной передаточной функции
num = [b];
den = [1, a1, a2];
Wz = tf(num, den, Td);  % Ts - интервал дискретизации

% Расчет отклика
y1 = lsim(Wz, u1, t);
y2 = lsim(Wz, u2, t);

gamma = 1;

figure;

% --- идентификация a,b ---
theta1 = zeros(3, N);      % [a1_hat; a2_hat; b_hat]
theta_hat1 = [0; 0; 0];       % начальное значение
phi1 = zeros(3,2);         % [y(k); u(k)]

for k_i = 3:N
    phi1 = [-y1(k_i-1); -y1(k_i-2); u1(k_i-2)];

    % e0(k)
    e0 = y1(k_i) - phi1.'*theta_hat1;

    % обновление
    theta_hat1 = theta_hat1 + gamma * (phi1 * e0)/(1 + gamma*(phi1.'*phi1));

    theta1(:, k_i) = theta_hat1;
end

% --- идентификация a,b ---
theta2 = zeros(3, N);      % [a1_hat; a2_hat; b_hat]
theta_hat2 = [0; 0; 0];       % начальное значение
phi2 = zeros(3,2);         % [y(k); u(k)]

for k_i = 3:N
    phi2 = [-y2(k_i-1); -y2(k_i-2); u2(k_i-2)];

    % e0(k)
    e0 = y2(k_i) - phi2.'*theta_hat2;

    % обновление
    theta_hat2 = theta_hat2 + gamma * (phi2 * e0)/(1 + gamma*(phi2.'*phi2));

    theta2(:, k_i) = theta_hat2;
end

% ---- графики ----
subplot(3, 1, 1);
plot(t, theta1(1,:), 'LineWidth', 1.6); hold on;
plot(t, theta2(1,:), 'LineWidth', 1.6);
yline(a1, '--k', 'LineWidth', 1.5); grid on;
legend(["одна гармоника", "две гармоники"], 'Location', 'northeast')
title('Оценки a_1');
xlabel("t")
ylabel("\theta(t)")
yticks([a1 -1 0])
ylim([-2 0.1])

subplot(3, 1, 2);
plot(t, theta1(2,:), 'LineWidth', 1.6); hold on;
plot(t, theta2(2,:), 'LineWidth', 1.6); hold on;
yline(a2, '--k', 'LineWidth', 1.5); grid on;
title('Оценки a_2');
legend(["одна гармоника", "две гармоники"], 'Location', 'southeast')
yticks([0 0.5 a2])
xlabel("t")
ylabel("\theta(t)")
ylim([-0.1 1]);

subplot(3, 1, 3);
plot(t, theta1(3,:), 'LineWidth', 1.6); hold on;
plot(t, theta2(3,:), 'LineWidth', 1.6); hold on;
yline(b, '--k', 'LineWidth', 1.5); grid on;
title('Оценки b');
ylim([-0.1 2.7])
yticks([0 1 2 b])
xlabel("t")
ylabel("\theta(t)")
legend(["одна гармоника", "две гармоники"], 'Location', 'southeast')
saveas(gcf, "images/task2_4.png")

% ---- нормы ошибок ----
err1 = sqrt((theta1(1,:) - a1).^2 + (theta1(2,:) - a2).^2 + (theta1(3,:) - b).^2);
err2 = sqrt((theta2(1,:) - a1).^2 + (theta2(2,:) - a2).^2 + (theta2(3,:) - b).^2);

figure;
plot(t, err1, 'LineWidth', 1.6); hold on;
plot(t, err2, 'LineWidth', 1.6);
grid on;
yline(0, '--k', 'LineWidth', 1.5);
legend(["одна гармоника", "две гармоники"], 'Location', 'northeast');
title("Нормы параметрических ошибок");
ylabel("f(t)");
xlabel("t");
ylim([-0.1 3.5])
saveas(gcf, "images/task2_4_errors.png")
