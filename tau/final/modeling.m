%% Параметры системы
n = 11;
rng(n, "philox");
M = randi([100000 1000000]) / 1000 / sqrt(2);
m = randi([1000 10000]) / 1000 * sqrt(3);
l = randi([100 1000]) / sqrt(5) / 100;
g = 9.81;

fprintf('Параметры системы:\n');
fprintf('M = %.4f кг (масса тележки)\n', M);
fprintf('m = %.4f кг (масса маятника)\n', m);
fprintf('l = %.4f м (длина маятника)\n', l);
fprintf('g = %.4f м/с² (ускорение свободного падения)\n\n', g);

% Линеаризованные матрицы системы
A = [0, 1, 0, 0;
    0, 0, 3*m*g/(4*M + m), 0; 
    0, 0, 0, 1; 
    0, 0, 6*(M + m)*g/l/(4*M + m), 0];

B = [0; 
    4 / (4*M + m);
    0; 
    6/l/(4*M + m)];

C = [1, 0, 0, 0;
    0, 0, 1, 0];

D = [0; 
    6/l/(4*M + m); 
    0; 
    12*(M + m) / m / l^2 / (4*M + m)];

% Нелинейная модель системы "маятник на тележке"
function dxdt = nonlinear_pendulum(t, x, u, f, M, m, l, g)
    % x = [x1; x2; x3; x4] = [позиция тележки; скорость тележки; угол маятника; угловая скорость]
    % u - управляющее воздействие на тележку
    % f - внешняя сила на маятник
    
    x1 = x(1); % позиция тележки
    x2 = x(2); % скорость тележки  
    x3 = x(3); % угол маятника
    x4 = x(4); % угловая скорость маятника
    
    % Вычисляем управляющие воздействия (если они функции времени)
    u_val = u;
    f_val = f;
    
    % Нелинейные уравнения движения
    denominator1 = (1/3)*(M + m)*l - (1/4)*m*l*cos(x3)^2;
    numerator1 = (1/3)*l*(u_val - 0.5*m*l*x4^2*sin(x3)) + 0.5*cos(x3)*(f_val + 0.5*m*g*l*sin(x3));
    dx2 = numerator1 / denominator1;
    
    denominator2 = (1/3)*(M + m)*m*l^2 - (1/4)*m^2*l^2*cos(x3)^2;
    numerator2 = 0.5*m*l*cos(x3)*(u_val - 0.5*m*l*x4^2*sin(x3)) + (M + m)*(f_val + 0.5*m*g*l*sin(x3));
    dx4 = numerator2 / denominator2;
    
    dxdt = [x2; dx2; x4; dx4];
end

% Параметры моделирования
t_span = [0 4];  % временной интервал
dt = 0.001;        % шаг интегрирования
t_values = 0:dt:4;

% Начальные условия из CS_project (все состояния ненулевые)
x01 = [0.001; 0.001; 0; 0];
x02 = [0; 0; 0.001; 0.001];
x03 = [-0.001; -0.001; -0.001; -0.001];
x04 = [-0.005; -0.003; 0.007; -0.002];



%% Управляющие воздействия
u = 0;  % без управления
f = 0;  % без внешней силы

fprintf('Решение 1: малые положительные значения [0.001, 0.001, 0, 0]\n');
[t1, sol1] = ode45(@(t,x) nonlinear_pendulum(t, x, u, f, M, m, l, g), t_values, x01);

fprintf('Решение 2: отрицательные значения [0, 0, 0.001, 0.001]\n');
[t2, sol2] = ode45(@(t,x) nonlinear_pendulum(t, x, u, f, M, m, l, g), t_values, x02);

fprintf('Решение 3: очень малые значения [-0.001, -0.001, -0.001, -0.001]\n');
[t3, sol3] = ode45(@(t,x) nonlinear_pendulum(t, x, u, f, M, m, l, g), t_values, x03);

fprintf('Решение 4: отрицательные значения [-0.005, -0.003, 0.007, -0.002]\n');
[t4, sol4] = ode45(@(t,x) nonlinear_pendulum(t, x, u, f, M, m, l, g), t_values, x04);

fprintf('\n=== МОДЕЛИРОВАНИЕ ЛИНЕАРИЗОВАННОЙ СИСТЕМЫ ===\n\n');

% Линеаризованная система: dx/dt = A*x + B*u
% Для сравнения используем те же начальные условия
fprintf('Линейная модель - Решение 1: [0.001, 0.001, 0, 0]\n');
[t1_lin, sol1_lin] = ode45(@(t,x) A*x + B*u, t_values, x01);

fprintf('Линейная модель - Решение 2: [0, 0, 0.001, 0.001]\n');
[t2_lin, sol2_lin] = ode45(@(t,x) A*x + B*u, t_values, x02);

fprintf('Линейная модель - Решение 3: [-0.001, -0.001, -0.001, -0.001]\n');
[t3_lin, sol3_lin] = ode45(@(t,x) A*x + B*u, t_values, x03);

fprintf('Линейная модель - Решение 4: [-0.005, -0.003, 0.007, -0.002]\n');
[t4_lin, sol4_lin] = ode45(@(t,x) A*x + B*u, t_values, x04);

fprintf('\n=== ВИЗУАЛИЗАЦИЯ РЕШЕНИЙ ===\n\n');

% Названия компонент состояния
state_vars = {"x_1", "x_2", "x_3", "x_4"};


%%
% Решение 1: Сравнение нелинейной и линейной моделей
fprintf('Создание графиков для решения 1...\n');
for i = 1:4
    figure;
    plot(t1, sol1(:,i), 'LineWidth', 3, 'DisplayName', 'Нелинейная модель', "Clipping", "on");
    hold on;
    plot(t1_lin, sol1_lin(:,i), 'LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Линейная модель', "Clipping", "on");
    xlabel('t', 'FontSize', 15);
    ylabel("x", 'FontSize', 15);
    title(state_vars{i} + "(t) при большом времени", 'FontSize', 17);
    legend('Location', 'northwest', 'FontSize', 13);
    set(gca, 'FontSize', 14);
    grid on;
    
    % Добавляем 10% отступы сверху и снизу
    y_combined = [sol1(:,i); sol1_lin(:,i)];
    y_min = min(y_combined);
    y_max = max(y_combined);
    y_range = y_max - y_min;
    margin = max(0.1 * y_range, 1e-6);
    ylim([y_min - margin, y_max + margin]);

    
    % Сохранение графика
    filename = sprintf('images/solution1_%s.png', state_vars{i});
    saveas(gcf, filename);
    fprintf('Сохранен график: %s\n', filename);
    close(gcf);
end
%%
% Решение 2: Сравнение нелинейной и линейной моделей
fprintf('Создание графиков для решения 2...\n');
for i = 1:4
    figure;
    plot(t2, sol2(:,i), 'LineWidth', 3, 'DisplayName', 'Нелинейная модель', "Clipping", "on");
    hold on;
    plot(t2_lin, sol2_lin(:,i), 'LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Линейная модель', "Clipping", "on");
    xlabel('t', 'FontSize', 15);
    ylabel("x", 'FontSize', 15);
    title(state_vars{i} + "(t) при большом времени", 'FontSize', 17);
    legend('Location', 'northwest', 'FontSize', 13);
    set(gca, 'FontSize', 14);
    grid on;
    
    % Добавляем 10% отступы сверху и снизу
    y_combined = [sol2(:,i); sol2_lin(:,i)];
    y_min = min(y_combined);
    y_max = max(y_combined);
    y_range = y_max - y_min;
    margin = max(0.1 * y_range, 1e-6);
    ylim([y_min - margin, y_max + margin]);
    
    % Сохранение графика
    filename = sprintf('images/solution2_%s.png', state_vars{i});
    saveas(gcf, filename);
    fprintf('Сохранен график: %s\n', filename);
    close(gcf);
end

% Решение 3: Сравнение нелинейной и линейной моделей
fprintf('Создание графиков для решения 3...\n');
for i = 1:4
    figure;
    plot(t3, sol3(:,i), 'LineWidth', 3, 'DisplayName', 'Нелинейная модель', "Clipping", "on");
    hold on;
    plot(t3_lin, sol3_lin(:,i), 'LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Линейная модель', "Clipping", "on");
    xlabel('t', 'FontSize', 15);
    ylabel("x", 'FontSize', 15);
    title(state_vars{i} + "(t) при большом времени", 'FontSize', 17);
    legend('Location', 'southwest', 'FontSize', 13);
    set(gca, 'FontSize', 14);
    grid on;
    
    % Добавляем 10% отступы сверху и снизу
    y_combined = [sol3(:,i); sol3_lin(:,i)];
    y_min = min(y_combined);
    y_max = max(y_combined);
    y_range = y_max - y_min;
    margin = max(0.1 * y_range, 1e-6);
    ylim([y_min - margin, y_max + margin]);
    
    % Сохранение графика
    filename = sprintf('images/solution3_%s.png', state_vars{i});
    saveas(gcf, filename);
    fprintf('Сохранен график: %s\n', filename);
    close(gcf);
end

% Решение 4: Сравнение нелинейной и линейной моделей
fprintf('Создание графиков для решения 4...\n');
for i = 1:4
    figure;
    plot(t4, sol4(:,i), 'LineWidth', 3, 'DisplayName', 'Нелинейная модель', "Clipping", "on");
    hold on;
    plot(t4_lin, sol4_lin(:,i), 'LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Линейная модель', "Clipping", "on");
    xlabel('t', 'FontSize', 15);
    ylabel("x", 'FontSize', 15);
    title(state_vars{i} + "(t) при большом времени", 'FontSize', 17);
    legend('Location', 'northwest', 'FontSize', 13);
    set(gca, 'FontSize', 14);
    grid on;
    
    % Добавляем 10% отступы сверху и снизу
    y_combined = [sol4(:,i); sol4_lin(:,i)];
    y_min = min(y_combined);
    y_max = max(y_combined);
    y_range = y_max - y_min;
    margin = max(0.1 * y_range, 1e-6);
    ylim([y_min - margin, y_max + margin]);

    
    % Сохранение графика
    filename = sprintf('images/solution4_%s.png', state_vars{i});
    saveas(gcf, filename);
    fprintf('Сохранен график: %s\n', filename);
    close(gcf);
end

%%
% Дополнительные графики с ограничением времени до 1 секунд
fprintf('\n=== СОЗДАНИЕ ГРАФИКОВ С ВРЕМЕНЕМ ДО 1 СЕКУНД ===\n\n');

% Решение 1: Графики с xlim(0, 1)
fprintf('Создание графиков для решения 1 (0-1 сек)...\n');
for i = 1:4
    figure;
    plot(t1, sol1(:,i), 'LineWidth', 3, 'DisplayName', 'Нелинейная модель', "Clipping", "on");
    hold on;
    plot(t1_lin, sol1_lin(:,i), 'LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Линейная модель', "Clipping", "on");
    xlabel('t', 'FontSize', 15);
    ylabel("x", 'FontSize', 15);
    title(state_vars{i} + "(t) при малом времени", 'FontSize', 17);
    legend('Location', 'northwest', 'FontSize', 13);
    set(gca, 'FontSize', 14);
    grid on;
    xlim([0, 1]);

    
    % Сохранение графика
    filename = sprintf('images/solution1_%s_1sec.png', state_vars{i});
    saveas(gcf, filename);
    fprintf('Сохранен график: %s\n', filename);
    close(gcf);
end
%%
% Решение 2: Графики с xlim(0, 1)
fprintf('Создание графиков для решения 2 (0-1 сек)...\n');
i = 1;
figure;
plot(t2, sol2(:,i), 'LineWidth', 3, 'DisplayName', 'Нелинейная модель', "Clipping", "on");
hold on;
plot(t2_lin, sol2_lin(:,i), 'LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Линейная модель', "Clipping", "on");
xlabel('t', 'FontSize', 15);
ylabel("x", 'FontSize', 15);
title(state_vars{i} + "(t) при малом времени", 'FontSize', 17);
legend('Location', 'northwest', 'FontSize', 13);
set(gca, 'FontSize', 14);
grid on;
xlim([0, 1]);
ylim([-0.25*1e-5, 3.5*1e-5]);
% Сохранение графика
filename = sprintf('images/solution2_%s_1sec.png', state_vars{i});
saveas(gcf, filename);
fprintf('Сохранен график: %s\n', filename);
close(gcf);
i=2;
figure;
plot(t2, sol2(:,i), 'LineWidth', 3, 'DisplayName', 'Нелинейная модель', "Clipping", "on");
hold on;
plot(t2_lin, sol2_lin(:,i), 'LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Линейная модель', "Clipping", "on");
xlabel('t', 'FontSize', 15);
ylabel("x", 'FontSize', 15);
title(state_vars{i} + "(t) при малом времени", 'FontSize', 17);
legend('Location', 'northwest', 'FontSize', 13);
set(gca, 'FontSize', 14);
grid on;
xlim([0, 1]);
ylim([-0.05*1e-4, 1.2*1e-4]);
% Сохранение графика
filename = sprintf('images/solution2_%s_1sec.png', state_vars{i});
saveas(gcf, filename);
fprintf('Сохранен график: %s\n', filename);
close(gcf);
for i = 3:4
    figure;
    plot(t2, sol2(:,i), 'LineWidth', 3, 'DisplayName', 'Нелинейная модель', "Clipping", "on");
    hold on;
    plot(t2_lin, sol2_lin(:,i), 'LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Линейная модель', "Clipping", "on");
    xlabel('t', 'FontSize', 15);
    ylabel("x", 'FontSize', 15);
    title(state_vars{i} + "(t) при малом времени", 'FontSize', 17);
    legend('Location', 'northwest', 'FontSize', 13);
    set(gca, 'FontSize', 14);
    grid on;
    xlim([0, 1]);

    % Сохранение графика
    filename = sprintf('images/solution2_%s_1sec.png', state_vars{i});
    saveas(gcf, filename);
    fprintf('Сохранен график: %s\n', filename);
    close(gcf);
end

%%
% Решение 3: Графики с xlim(0, 1)
fprintf('Создание графиков для решения 3 (0-1 сек)...\n');
for i = 1:4
    figure;
    plot(t3, sol3(:,i), 'LineWidth', 3, 'DisplayName', 'Нелинейная модель', "Clipping", "on");
    hold on;
    plot(t3_lin, sol3_lin(:,i), 'LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Линейная модель', "Clipping", "on");
    xlabel('t', 'FontSize', 15);
    ylabel("x", 'FontSize', 15);
    title(state_vars{i} + "(t) при малом времени", 'FontSize', 17);
    legend('Location', 'southwest', 'FontSize', 13);
    set(gca, 'FontSize', 14);
    grid on;
    xlim([0, 1]);
    
    % Сохранение графика
    filename = sprintf('images/solution3_%s_1sec.png', state_vars{i});
    saveas(gcf, filename);
    fprintf('Сохранен график: %s\n', filename);
    close(gcf);
end

% Решение 4: Графики с xlim(0, 1)
fprintf('Создание графиков для решения 4 (0-1 сек)...\n');
i = 1;
figure;
plot(t4, sol4(:,i), 'LineWidth', 3, 'DisplayName', 'Нелинейная модель', "Clipping", "on");
hold on;
plot(t4_lin, sol4_lin(:,i), 'LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Линейная модель', "Clipping", "on");
xlabel('t', 'FontSize', 15);
ylabel("x", 'FontSize', 15);
title(state_vars{i} + "(t) при малом времени", 'FontSize', 17);
legend('Location', 'northeast', 'FontSize', 13);
set(gca, 'FontSize', 14);
grid on;
xlim([0, 1]);
% Сохранение графика
filename = sprintf('images/solution4_%s_1sec.png', state_vars{i});
saveas(gcf, filename);
fprintf('Сохранен график: %s\n', filename);
close(gcf);
for i = 2:4
    figure;
    plot(t4, sol4(:,i), 'LineWidth', 3, 'DisplayName', 'Нелинейная модель', "Clipping", "on");
    hold on;
    plot(t4_lin, sol4_lin(:,i), 'LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Линейная модель', "Clipping", "on");
    xlabel('t', 'FontSize', 15);
    ylabel("x", 'FontSize', 15);
    title(state_vars{i} + "(t) при малом времени", 'FontSize', 17);
    legend('Location', 'northwest', 'FontSize', 13);
    set(gca, 'FontSize', 14);
    grid on;
    xlim([0, 1]);

    % Сохранение графика
    filename = sprintf('images/solution4_%s_1sec.png', state_vars{i});
    saveas(gcf, filename);
    fprintf('Сохранен график: %s\n', filename);
    close(gcf);
end

fprintf('\nВсе графики сохранены в папку images/\n');
fprintf('Создано 24 отдельных графика:\n');
fprintf('- 12 графиков для полного времени (4 сек)\n');
fprintf('- 12 графиков для начального периода (1 сек)\n');

