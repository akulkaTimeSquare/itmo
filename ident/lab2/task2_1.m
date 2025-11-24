load('ident_lab2_v05.mat');
% исходные параметры
a1 = zad2.a1;
a2 = zad2.a2;
b = zad2.b;
w = zad2.w;

% Параметры
Td = 0.1;
t_end = 60;
t = 0:Td:t_end;
N = t_end/Td+1;

% Входной сигнал
u = sin(w*t);

% Создание дискретной передаточной функции
num = [b];
den = [1, a1, a2];
Wz = tf(num, den, Td);  % Ts - интервал дискретизации

% Расчет отклика
y = lsim(Wz, u, t);

gamma = 1;

figure;

% --- идентификация a,b ---
theta = zeros(3, N);      % [a1_hat; a2_hat; b_hat]
theta_hat = [0; 0; 0];       % начальное значение
phi = zeros(3,2);         % [y(k); u(k)]

for k_i = 3:N
    phi = [-y(k_i-1); -y(k_i-2); u(k_i-2)];

    % e0(k)
    e0 = y(k_i) - phi.'*theta_hat;

    % обновление
    theta_hat = theta_hat + gamma * (phi * e0)/(1 + gamma*(phi.'*phi));

    theta(:, k_i) = theta_hat;
end

% ---- графики ----
subplot(3, 1, 1);
plot(t, theta(1,:), 'LineWidth', 1.6); hold on;
yline(a1, '--k', 'LineWidth', 1.5); grid on;
title('Оценка a_1 при u = sin(\omegat)');
yticks([a1 -1 0])
xlabel("t")
ylabel("\theta(t)")
ylim([-2 0.1])

subplot(3, 1, 2);
plot(t, theta(2,:), 'LineWidth', 1.6); hold on;
yline(a2, '--k', 'LineWidth', 1.5); grid on;
title('Оценка a_2 при u = sin(\omegat)');
yticks([0 0.5 a2])
xlabel("t")
ylabel("\theta(t)")
ylim([-0.1 1]);

subplot(3, 1, 3);
plot(t, theta(3,:), 'LineWidth', 1.6); hold on;
yline(b, '--k', 'LineWidth', 1.5); grid on;
title('Оценка b при u = sin(\omegat)');
ylim([-0.2 2.7])
yticks([0 1 2 b])
xlabel("t")
ylabel("\theta(t)")
saveas(gcf, "images/task2_1.png")