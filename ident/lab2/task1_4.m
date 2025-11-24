load('ident_lab2_v05.mat');
% ---- параметры системы ----
a = zad1.a;
b = zad1.b;
w = zad1.w;

Td = 0.1;
Tend = 30;

k = 0:Td:Tend;
N = length(k);

% вход
u = sin(w * k);

% Создание дискретной передаточной функции
num = [b];
den = [1, a];
Wz = tf(num, den, Td);  % Ts - интервал дискретизации

% Расчет отклика
y = lsim(Wz, u, k);

% ---- алгоритмы для gamma ----
gammas = [0.5, 10];

figure;

gi = 1;
gamma = gammas(gi);

% параметры оценки
theta = zeros(2, N);      % [a_hat; b_hat]
theta_hat = [0; 0];       % начальная оценка

for k_i = 2:N
    phi = [-y(k_i-1); u(k_i-1)];

    % ошибка
    e0 = y(k_i) - phi.'*theta_hat;

    % градиентный алгоритм
    theta_hat = theta_hat + gamma * phi * e0;

    theta(:, k_i) = theta_hat;
end

subplot(2,1,1);
plot(k, theta(1,:), 'LineWidth', 1.5); hold on;
yline(a, '--k', 'LineWidth', 1.5); grid on;
yticks([-0.8 -0.5, 0, 0.5, a]);
ylim([-0.8 1.1])
xlabel("t")
ylabel("\theta(t)")
title(['Оценка a при \gamma = ', num2str(gamma)]);

subplot(2,1,2);
plot(k, theta(2,:), 'LineWidth', 1.5); hold on;
yline(b, '--k', 'LineWidth', 1.5); grid on;
yticks([0 1 2 b])
title(['Оценка b при \gamma = ', num2str(gamma)]);
set(gcf, 'PaperPositionMode','auto');
xlabel("t")
ylabel("\theta(t)")
saveas(gcf, "images/task1_4.png")
