load('ident_lab2_v05.mat');
a = zad3.a;
b = zad3.b;
w = zad3.w;
theta_true = [a; b];

dt = 0.0001;
Tfinal = 15;
t = 0:dt:Tfinal;
N = length(t);

u = sin(w*t);
y = zeros(1,N);
y(1) = 0;

% --- estimator parameters ---
gammas = [1, 3, 10];
i = 2;
gamma = gammas(i);
theta_hat = zeros(2, N);   % theta_hat(:,k,i) — оценка в момент k для i-го gamma

% --- симуляция (Euler для системы и для адаптации) ---
for k = 1:N-1
    % --- истинная динамика: вычислим производную (доступна измерению) ---
    y_dot = -a * y(k) + b * u(k);   % т.к. у нас идеальная модель
    % интегрируем plant (Euler)
    y(k+1) = y(k) + dt * y_dot;
    
    % регрессор
    phi = [-y(k); u(k)];              % (2x1)
    
    th = theta_hat(:,k);
    ydot_pred = phi' * th;        % предсказание производной
    e = y_dot - ydot_pred;        % ошибка (скаляр)
    % адаптационный закон (интегрируем тоже Euler'ом)
    th_dot = gamma * phi * e; % (2x1)
    theta_hat(:,k+1) = th + dt * th_dot;
end

% --- Графики: оценки параметров против истинных значений ---
figure;
% a (theta(1))
subplot(2,1,1);
hold on; grid on;
colors = lines(nG);
plot(t, theta_hat(1,:), 'LineWidth', 1.5);
yline(a, '--k', 'LineWidth', 1.5); % истинное a0
title(['Оценка a при \gamma = ', num2str(gamma)]);
xlabel('t'); ylabel('\theta(t)');
ylim([-0.7 0.6])
yticks([-0.7 -0.35 0 0.35 0.5])
hold off;

% b (theta(2))
subplot(2,1,2);
hold on; grid on;
plot(t, theta_hat(2,:), 'LineWidth', 1.5);
yline(b, '--k', 'LineWidth',1.5); % истинное b
title(['Оценка b при \gamma = ', num2str(gamma)]);
set(gcf, 'PaperPositionMode','auto');
xlabel('t'); ylabel('\theta(t)');
hold off;
yticks([0 1 2 b 3.5])
ylim([0 3.5])
saveas(gcf, "images/task3_2.png")