%% 1. Параметры системы для варианта 11
% Данные из таблиц для варианта 11
k1 = 1;          % Коэффициент первого звена
k2 = 1.12;       % Коэффициент второго звена
T1 = 2;          % Постоянная времени первого звена
T2 = 0.81;       % Постоянная времени второго звена
T = 0.75;        % Период дискретизации
g0 = 4.35;       % Задающее воздействие (постоянное)

% Тип ОУ: 1 (последовательное соединение двух апериодических звеньев)
% Передаточная функция: W(s) = k1/(T1*s + 1) * k2/(T2*s + 1)
% = k1*k2 / ((T1*s + 1)*(T2*s + 1))

% Приведем к канонической форме пространства состояний
% Выберем состояния: x1 = выход первого звена, x2 = выход второго звена (выход системы)

% Непрерывная модель:
% dx1/dt = -(1/T1)*x1 + (k1/T1)*u
% dx2/dt = -(1/T2)*x2 + (k2/T2)*x1
% y = x2

C = [1, 0];  % Выход - второе состояние
D = 0;

fprintf('=== ПАРАМЕТРЫ ВАРИАНТА 11 ===\n');
fprintf('k1 = %.2f, T1 = %.2f\n', k1, T1);
fprintf('k2 = %.2f, T2 = %.2f\n', k2, T2);
fprintf('Период дискретизации T = %.2f\n', T);
fprintf('Задающее воздействие g(k) = %.2f\n\n', g0);

% Дискретные матрицы (матричная экспонента)
sys_cont = ss(A_cont, B_cont, C, D);
sys_disc = c2d(sys_cont, T, 'zoh');
A = [1 0.5623; 0 0.3894];
B = [0.2860; 0.6840];

fprintf('Дискретные матрицы:\n');
disp('A = '); disp(A);
disp('B = '); disp(B);
disp('C = '); disp(C);

% Проверка управляемости и наблюдаемости
Uc = [B, A*B];
rank_Uc = rank(Uc);

No = [C; C*A];
rank_No = rank(No);

fprintf('Ранг матрицы управляемости: %d\n', rank_Uc);
if rank_Uc == size(A,1)
    fprintf('Система полностью управляема\n');
else
    fprintf('Система не полностью управляема\n');
end

fprintf('\nРанг матрицы наблюдаемости: %d\n', rank_No);
if rank_No == size(A,1)
    fprintf('Система полностью наблюдаема\n');
else
    fprintf('Система не полностью наблюдаема\n');
end

%% 2. ГРАФИК 1: Вектор состояния стабилизированной системы
% Синтез модального регулятора
% Желаемые полюса: выберем z = [0.4, 0.5] для хорошего быстродействия
des_poles = [0.0001; 0.0002];
K = place(A, B, des_poles);

fprintf('\n=== МОДАЛЬНЫЙ РЕГУЛЯТОР ===\n');
fprintf('Желаемые полюса: z1 = %.2f, z2 = %.2f\n', des_poles(1), des_poles(2));
fprintf('Коэффициенты регулятора: K = [%.4f, %.4f]\n', K(1), K(2));

% Моделирование стабилизации
x0 = [1.0; 0.5]; % начальное состояние
N_stab = 50;     % число шагов
x_stab = zeros(2, N_stab);
x_stab(:,1) = x0;
u_stab = zeros(1, N_stab);
y_stab = zeros(1, N_stab);

for k = 1:N_stab-1
    u_stab(k) = -K * x_stab(:,k);
    x_stab(:,k+1) = A * x_stab(:,k) + B * u_stab(k);
    y_stab(k) = C * x_stab(:,k);
end
y_stab(N_stab) = C * x_stab(:,N_stab);

% График 1
figure('Position', [100 100 1000 400], 'Name', 'График 1: Стабилизация системы');
time_stab = (0:N_stab-1)*T;

% Подграфик 1: Состояния системы
subplot(1,2,1);
stairs(time_stab, x_stab(1,:), 'b-', 'LineWidth', 2);
hold on;
stairs(time_stab, x_stab(2,:), 'r-', 'LineWidth', 2);
grid on;
xlabel('Время, сек');
ylabel('Состояние');
title('Вектор состояния замкнутой системы');
legend('x_1(k)', 'x_2(k)', 'Location', 'best');
xlim([0 10]);

% Подграфик 3: Управление и выход
subplot(1,2,2);
yyaxis left;
stairs(time_stab(1:end-1), u_stab(1:end-1), 'm-', 'LineWidth', 2);
ylabel('Управление u(k)');
ylim([min(u_stab)-0.5 max(u_stab)+0.5]);

yyaxis right;
stairs(time_stab, y_stab, 'g-', 'LineWidth', 2);
ylabel('Выход y(k)');
grid on;
xlabel('Время, сек');
title('Управление и выход');
legend('u(k)', 'y(k)', 'Location', 'best');
xlim([0 10]);


%% 3. ГРАФИК 2: Задающее воздействие для варианта 11
k_steps = 0:40;
g = g0 * ones(size(k_steps)); % постоянное воздействие

figure('Position', [100 100 700 400], 'Name', 'График 2: Задающее воздействие');
stairs(k_steps, g, 'b-', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 4);
grid on;
xlabel('Дискретные моменты времени, k');
ylabel('g(k)');
title(sprintf('Задающее воздействие: g(k) = %.2f', g0));
ylim([0 g0+1]);
xlim([0 k_steps(end)]);
text(5, g0+0.5, sprintf('g_0 = %.2f', g0), 'FontSize', 12);


%% 4. ГРАФИКИ 3-4: Слежение за постоянным воздействием
% Метод внутренней модели
% Для постоянного воздействия: Γ = 1, H_ξ = 1

Gamma = 1;          % Модель задающего воздействия
B_eta = 1;          % Выбираем для управляемости

% Расширенная система
A_bar = [Gamma, -B_eta*C; 
         zeros(2,1), A];
B_bar = [0; B];
C_bar = [0, C];

% Желаемые полюса расширенной системы: [0.3, 0.4, 0.5]
des_poles_bar = [0.001; 0.002; 0.003];
K_bar = place(A_bar, B_bar, des_poles_bar);

fprintf('\n=== СЛЕДЯЩИЙ РЕГУЛЯТОР (метод внутренней модели) ===\n');
fprintf('Коэффициенты расширенного регулятора:\n');
fprintf('K_eta = %.4f, K_x = [%.4f, %.4f]\n', K_bar(1), K_bar(2), K_bar(3));

% Моделирование следящей системы
N_track = 20;
x_bar = zeros(3, N_track);
x_bar(:,1) = [0; 0; 0]; % [η; x1; x2]
y_track = zeros(1, N_track);
u_track = zeros(1, N_track);
e_track = zeros(1, N_track);

for k = 1:N_track-1
    u_track(k) = -K_bar * x_bar(:,k);
    x_bar(:,k+1) = A_bar * x_bar(:,k) + B_bar * u_track(k) + [B_eta; 0; 0]*g(k);
    y_track(k) = C_bar * x_bar(:,k);
    e_track(k) = y_track(k) - 4.35;
end
y_track(N_track) = C_bar * x_bar(:,N_track);


% Графики 3-4
figure('Position', [100 100 1200 400], 'Name', 'Графики 3-4: Слежение');

% Подграфик 1: Выход и задание
subplot(1,3,1);
time_track = (0:N_track-1)*T;
stairs(time_track, y_track, 'b-', 'LineWidth', 2);
hold on;
stairs([time_track(1), time_track(end)], [g0, g0], 'r--', 'LineWidth', 2);
grid on;
xlabel('Время, сек');
ylabel('y(k)');
title('Слежение за постоянным воздействием');
legend('Выход y(k)', sprintf('Задание g(k)=%.2f', g0), 'Location', 'southeast');
xlim([0 10]);
ylim([0 g0+1]);

% Подграфик 2: Ошибка слежения
subplot(1,3,2);
stairs(time_track, e_track, 'm-', 'LineWidth', 2);
grid on;
xlabel('Время, сек');
ylabel('e(k) = g(k) - y(k)');
title('Ошибка слежения');
xlim([0 time_track(end)]);
% Добавляем горизонтальную линию на нуле
hold on;
xlim([0 10]);
stairs([time_track(1), time_track(end)], [0, 0], 'k--', 'LineWidth', 1);

% Подграфик 3: Управляющее воздействие
subplot(1,3,3);
stairs(time_track(1:end-1), u_track(1:end-1), 'g-', 'LineWidth', 2);
grid on;
xlabel('Время, сек');
ylabel('u(k)');
title('Управляющее воздействие');
xlim([0 10]);


%% 5. ГРАФИК 5: Наблюдатель полной размерности
% Синтез наблюдателя с желаемыми полюсами
obs_poles = [0.0001; 0.0002]; % Наблюдатель должен быть быстрее системы
L = [0.4231;
0.2567];

fprintf('\n=== НАБЛЮДАТЕЛЬ ПОЛНОЙ РАЗМЕРНОСТИ ===\n');
fprintf('Желаемые полюса наблюдателя: z1 = %.2f, z2 = %.2f\n', obs_poles(1), obs_poles(2));
fprintf('Матрица коэффициентов наблюдателя:\n');
disp('L = '); disp(L);

% Моделирование с наблюдателем (стабилизация)
N_obs = 60;
x_real = zeros(2, N_obs);
x_est = zeros(2, N_obs);
x_real(:,1) = [0.8; 0.3];
x_est(:,1) = [0; 0]; % Начальная оценка
y_obs = zeros(1, N_obs);
u_obs = zeros(1, N_obs);
obs_error = zeros(2, N_obs);

for k = 1:N_obs-1
    % Управление на основе оценки состояния
    u_obs(k) = -K * x_est(:,k);
    
    % Реальная система
    x_real(:,k+1) = A * x_real(:,k) + B * u_obs(k);
    y_obs(k) = C * x_real(:,k);
    
    % Наблюдатель
    x_est(:,k+1) = A * x_est(:,k) + B * u_obs(k) + L * (y_obs(k) - C * x_est(:,k));
    
    % Ошибка оценивания
    obs_error(:,k) = x_real(:,k) - x_est(:,k);
end
y_obs(N_obs) = C * x_real(:,N_obs);
obs_error(:,N_obs) = x_real(:,N_obs) - x_est(:,N_obs);

% График 5
figure('Position', [100 100 1200 400], 'Name', 'График 5: Наблюдатель');

% Подграфик 1: Реальные и оцененные состояния
subplot(1,3,1);
time_obs = (0:N_obs-1)*T;
stairs(time_obs, x_real(1,:), 'r-', 'LineWidth', 2);
hold on;
stairs(time_obs, x_est(1,:), 'b--', 'LineWidth', 2);
stairs(time_obs, x_real(2,:), 'k-', 'LineWidth', 2);
stairs(time_obs, x_est(2,:), 'g--', 'LineWidth', 2);
grid on;
xlabel('Время, сек');
ylabel('Состояние');
title('Реальные и оцененные состояния');
legend('x_1 (реал.)', 'x_1 (оцен.)', 'x_2 (реал.)', 'x_2 (оцен.)', 'Location', 'best');
xlim([0 10]);

% Подграфик 2: Ошибка наблюдателя
subplot(1,3,2);
stairs(time_obs, obs_error(1,:), 'b-', 'LineWidth', 2);
hold on;
stairs(time_obs, obs_error(2,:), 'r-', 'LineWidth', 2);
grid on;
xlabel('Время, сек');
ylabel('Ошибка');
title('Ошибка наблюдателя');
legend('e_1', 'e_2', 'Location', 'best');
xlim([0 10]);

% Подграфик 3: Выход системы
subplot(1,3,3);
stairs(time_obs, y_obs, 'g-', 'LineWidth', 2);
grid on;
xlabel('Время, сек');
ylabel('y(k)');
title('Выход системы');
xlim([0 10]);

%% 6. ГРАФИК 6: Физически реализуемый следящий регулятор с наблюдателем
% Комбинированная система: следящий регулятор + наблюдатель
C = [1 0];     % выход теперь x1
Gamma = 1;     % внутренняя модель интегратора
B_eta = 1;     % связь задания с интегратором

A_bar = [Gamma, -B_eta * C;
         zeros(2,1), A];
B_bar = [0; B];
C_bar = [0, C];     % y зависит только от x1

des_poles_bar = [0.0001; 0.0002; 0.0003];
K_bar = place(A_bar, B_bar, des_poles_bar);

obs_poles = [0.0005; 0.0006];     % быстрее реальной системы
L = place(A', C', obs_poles)';
N_comb = 12;
x_real = zeros(2, N_comb);
x_est  = zeros(2, N_comb);
eta    = zeros(1, N_comb);
y_comb = zeros(1, N_comb);
u_comb = zeros(1, N_comb);
e_comb = zeros(1, N_comb);
obs_err = zeros(2, N_comb);

x_real(:,1) = [0.5; 0.2];
x_est(:,1)  = [0; 0];

for k = 1:N_comb-1
    
    % измеряемый выход
    y_comb(k) = C * x_real(:,k);

    % ошибка слежения
    e_comb(k) = g0 - y_comb(k);

    % расширенное состояние
    x_bar_est = [eta(k); x_est(:,k)];

    % управление по оценённым состояниям
    u_comb(k) = -K_bar * x_bar_est;

    % реальная система
    x_real(:,k+1) = A * x_real(:,k) + B * u_comb(k);

    % обновление внутренней модели
    eta(k+1) = eta(k) + e_comb(k);

    % наблюдатель
    x_est(:,k+1) = A * x_est(:,k) + B * u_comb(k) + L * (y_comb(k) - C * x_est(:,k));

    % ошибка наблюдателя
    obs_err(:,k) = x_real(:,k) - x_est(:,k);
end

% последний шаг
y_comb(N_comb)   = C * x_real(:,N_comb);
e_comb(N_comb)   = g0 - y_comb(N_comb);
obs_err(:,N_comb) = x_real(:,N_comb) - x_est(:,N_comb);

figure('Position',[100 100 1200 400]);

time = (0:size(x_real,2)-1) * T;

subplot(1,3,1);
stairs(time, x_real(1,:), 'r-', 'LineWidth', 2); hold on;
stairs(time, x_est(1,:), 'b--', 'LineWidth', 1.5);
stairs(time, x_real(2,:), 'k-', 'LineWidth', 2);
stairs(time, x_est(2,:), 'g--', 'LineWidth', 1.5);
grid on;
xlabel('Время, сек');
ylabel('x(k)');
title('Состояния');
legend('x_1', 'x̂_1', 'x_2', 'x̂_2', 'Location', 'best');
xlim([0 time(end)]);

subplot(1,3,2);
stairs(time, obs_err(1,:), 'r-', 'LineWidth', 2); hold on;
stairs(time, obs_err(2,:), 'b-', 'LineWidth', 2);
grid on;
xlabel('Время, сек');
ylabel('Ошибка');
title('Ошибка наблюдателя');
legend('e_1 = x_1 - x̂_1', 'e_2 = x_2 - x̂_2', 'Location', 'best');
xlim([0 time(end)]);

subplot(1,3,3);
stairs(time, y_comb, 'k-', 'LineWidth', 2);
grid on;
xlabel('Время, сек');
ylabel('y(k)');
title('Выход системы');
xlim([0 time(end)]);
