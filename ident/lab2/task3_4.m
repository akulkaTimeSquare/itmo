load('ident_lab2_v05.mat');
a = zad3.a;
b = zad3.b;
w = zad3.w;
theta_true = [a; b];

dt = 0.0001;
Tfinal = 25;
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

% --- Графики: оценки параметров против истинных значений ---
figure;

% a (theta(1))
subplot(2,1,1);
hold on; grid on;
colors = lines(nG);
for i = 1:nG
    plot(t, squeeze(theta_hat(1,:,i)), 'LineWidth', 1.5, 'Color', colors(i,:));
end
yline(a, '--k', 'LineWidth',1.5); % истинное a0
title('Оценки a');
xlabel('t'); ylabel('\theta(t)');
ylim([-0.8 1.1])
yticks([-0.75 -0.5 -0.25 0 0.25 0.5 0.75 1])
legend([arrayfun(@(g) sprintf('\\gamma = %d',g), gammas, 'UniformOutput',false)], 'Location','southeast');
hold off;

% b (theta(2))
subplot(2,1,2);
hold on; grid on;
for i = 1:nG
    plot(t, squeeze(theta_hat(2,:,i)), 'LineWidth', 1.4, 'Color', colors(i,:));
end
yline(b, '--k', 'LineWidth',1.5); % истинное b
title('Оценки b');
xlabel('t, s'); ylabel('\theta(t)');
yticks([0 1 2 b 3.5])
ylim([0 3.5])
legend([arrayfun(@(g) sprintf('\\gamma = %d',g), gammas, 'UniformOutput',false)], 'Location','southeast');
set(gcf, 'PaperPositionMode','auto');
hold off;
saveas(gcf, "images/task3_4.png")

% --- График нормы ошибки оценок ---
figure;
hold on; grid on;
for i = 1:nG
    plot(t, err_norm(:,i), 'LineWidth', 1.5, 'Color', colors(i,:));
end
title('Нормы параметрических ошибок');
xlabel('t'); ylabel('f(t)');
yline(0, '--k', 'LineWidth', 1.5);
ylim([-0.1 3.5])
legend(arrayfun(@(g) sprintf('\\gamma = %d',g), gammas, 'UniformOutput',false), 'Location','northeast');
set(gcf, 'PaperPositionMode','auto');
hold off;
saveas(gcf, "images/task3_errors.png")