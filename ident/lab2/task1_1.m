load('ident_lab2_v05.mat');
% исходные параметры
a = zad1.a;
b = zad1.b;
w = zad1.w;

% Параметры
Td = 0.1;
t_end = 20;
t = 0:Td:t_end;
N = t_end/Td+1;

% Входной сигнал
u = sin(w*t);

% Создание дискретной передаточной функции
num = [b];
den = [1, a];
Wz = tf(num, den, Td);  % Ts - интервал дискретизации

% Расчет отклика
y = lsim(Wz, u, t);

gammas = [1, 3, 10];

figure;

gi = 1
gamma = gammas(gi);

theta = zeros(2, N);
theta_hat = [0; 0];
phi = zeros(2,1);
for k_i = 2:N
    phi = [-y(k_i-1); u(k_i-1)];
    e0 = y(k_i) - phi.'*theta_hat;
    theta_hat = theta_hat + gamma * (phi * e0)/(1 + gamma*(phi.'*phi));
    theta(:, k_i) = theta_hat;
end

subplot(2,1,1);
plot(t, theta(1,:), 'LineWidth', 1.5); hold on;
yline(a, '--k', 'LineWidth', 1.5); grid on;
yticks([-0.5, 0, 0.5, a]);
ylim([-0.5 1.1])
xlabel("t")
ylabel("\theta(t)")
title(['Оценка a при \gamma = ', num2str(gamma)]);

subplot(2,1,2);
plot(t, theta(2,:), 'LineWidth', 1.5); hold on;
yline(b, '--k', 'LineWidth', 1.5); grid on;
yticks([0 1 2 b])
title(['Оценка b при \gamma = ', num2str(gamma)]);
set(gcf, 'PaperPositionMode','auto');
xlabel("t")
ylabel("\theta(t)")
saveas(gcf, "images/task1_1.png")