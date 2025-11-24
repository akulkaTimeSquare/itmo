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
nG = length(gammas);
theta_hat = zeros(2, N, nG);   % theta_hat(:,k,i) — оценка в момент k для i-го gamma

% --- запомним ошибки и нормы ---
err_norm = zeros(N, nG);

% --- симуляция (Euler для системы и для адаптации) ---
for k = 1:N-1
    % --- истинная динамика: вычислим производную (доступна измерению) ---
    y_dot = -a * y(k) + b * u(k);   % т.к. у нас идеальная модель
    % интегрируем plant (Euler)
    y(k+1) = y(k) + dt * y_dot;
    
    % регрессор
    phi = [-y(k); u(k)];              % (2x1)
    
    for i = 1:nG
        th = theta_hat(:,k,i);
        ydot_pred = phi' * th;        % предсказание производной
        e = y_dot - ydot_pred;        % ошибка (скаляр)
        % адаптационный закон (интегрируем тоже Euler'ом)
        th_dot = gammas(i) * phi * e; % (2x1)
        theta_hat(:,k+1,i) = th + dt * th_dot;
        
        % запомним норму ошибки оценки
        err_norm(k,i) = norm(theta_hat(:,k,i) - theta_true);
    end
end
% финальные нормы
for i=1:nG, err_norm(N,i) = norm(theta_hat(:,N,i) - theta_true); end

% Графики
figure;
stairs(t, u, 'LineWidth', 1.5);
hold on;
stairs(t, y, 'LineWidth', 1.5);
hold off;
grid on;
legend('Вход u(t)', 'Выход y(t)');
title('Отклик непрерывной системы');
xlabel('t');
ylabel('f(t)');
ylim([-1.1 1.5])
set(gcf, 'PaperPositionMode','auto');
saveas(gcf, "images/task3.png")
