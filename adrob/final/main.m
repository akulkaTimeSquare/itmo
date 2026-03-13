function main()
    clear; clc; close all;

    %% 1. Определение системы
    A = [0, 1; 9, 2];
    b = [9; 5];
    d = b;
    C = [3, 3];
    
    % Проверка управляемости и наблюдаемости
    disp('Rank U:'); disp(rank([b, A*b]));
    disp('Rank W:'); disp(rank([C; C*A]));

    %% 2. Синтез регулятора (LQR)
    % Находим K, чтобы сделать матрицу A - bK Гурвицевой
    Q_lqr = eye(2);
    r_lqr = 3;
    % Используем стандартную функцию lqr (аналог icare)
    [K, P_lqr, ~] = lqr(A, b, Q_lqr, r_lqr);
    
    Acl = A - b*K;
    disp('Eigenvalues of Closed Loop (A - bK):');
    disp(eig(Acl));

    %% 3. Параметры адаптации и возмущения
    % Из отчета: f(t) = 0.5*sin(9t + 2) + 0.5*sin(7t - 16)
    % Базис регрессора: [cos(9t), sin(9t), cos(7t), sin(7t)]
    omega_f = [9, 7]; 
    
    % Вычисление истинных параметров theta_true
    % 0.5*sin(9t+2) = 0.5(sin(9t)cos(2) + cos(9t)sin(2))
    % Coeffs: cos(9t)->0.5*sin(2), sin(9t)->0.5*cos(2)
    % Аналогично для второй гармоники
    theta_true = [
        0.5 * sin(2);   % для cos(9t)
        0.5 * cos(2);   % для sin(9t)
        0.5 * sin(-16); % для cos(7t)
        0.5 * cos(-16)  % для sin(7t)
    ];

    %% 4. Настройка схемы Лайона
    num_harmonics = length(omega_f);
    m = 2 * num_harmonics; % Размерность вектора тета (4)
    n = size(A, 1);        % Размерность системы (2)

    lambda = [5, 10, 15, 20]; % Параметры фильтров
    q = length(lambda);
    gamma = 2000;            % Коэффициент скорости адаптации

    % Матрица P для уравнения Ляпунова (для наблюдателя/адаптации)
    Q_lyap = eye(2);
    P = lyap(Acl', Q_lyap);

    %% 5. Начальные условия
    x0 = [0; 0];
    theta_hat0 = zeros(m, 1);       
    z_filt0 = zeros(n, m);          % n x m
    z_theta_filt0 = zeros(n, 1);    % n x 1
    z_omega0 = zeros(q * m, 1);     % q*m x 1
    z_omega_theta0 = zeros(q, 1);   % q x 1
    z_e0 = zeros(q, 1);             % q x 1
    
    % Собираем полный вектор состояний
    % Порядок: x(2), theta_hat(4), z_filt(8), z_theta_filt(2), z_omega(16), z_omega_theta(4), z_e(4)
    y0_lyon = [x0; theta_hat0; z_filt0(:); z_theta_filt0; z_omega0; z_omega_theta0; z_e0];

    %% 6. Моделирование
    tspan = [0 10];
    opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-4);

    [t_lyon, y_lyon] = ode45(@(t,y) sim_lyon(t, y, A, b, K, d, P, omega_f, lambda, gamma, n, m, q, theta_true), ...
        tspan, y0_lyon, opts);

    %% 7. Извлечение и обработка результатов
    % Индексы для распаковки вектора состояний
    idx_x = 1:n;
    idx_theta = n+1 : n+m;
    
    x_res = y_lyon(:, idx_x);
    theta_hat_res = y_lyon(:, idx_theta);
    
    % Восстановление сигналов u(t), f(t), y(t)
    u_res = zeros(length(t_lyon), 1);
    f_res = zeros(length(t_lyon), 1);
    f_hat_res = zeros(length(t_lyon), 1);
    
    for i = 1:length(t_lyon)
        omega = generate_regressor(t_lyon(i), omega_f, m);
        
        % Управление: u = -Kx - theta_hat'*omega
        u_res(i) = -K * x_res(i, :)' - theta_hat_res(i, :) * omega;
        
        % Истинное возмущение
        f_res(i) = theta_true' * omega;
        
        % Оценка возмущения
        f_hat_res(i) = theta_hat_res(i, :) * omega;
    end
    
    y_out = (C * x_res')';
    error_f = f_res - f_hat_res;

    %% 8. Построение графиков
        % 1. Выход y(t)
    figure;
    plot(t_lyon, y_out, 'LineWidth', 2);
    title('Выходная переменная y(t) при модифицированном АА');
    xlabel('t'); ylabel('y(t)');
    ylim([-1.6 0.6]);
    grid on;
    saveas(gcf, 'images/lyon_y.png');

    % 2. Состояния x(t)
    figure;
    plot(t_lyon, x_res, 'LineWidth', 2);
    legend('x_1', 'x_2');
    title('Состояния системы x(t) при модифицированном АА');
    xlabel('t'); ylabel('x(t)');
    grid on;
    saveas(gcf, 'images/lyon_x.png');

    % 3. Параметры Theta
    figure;
    hold on;
    colors = lines(m);
    for i = 1:m
        plot(t_lyon, theta_hat_res(:, i), 'Color', colors(i,:), 'LineWidth', 2.0);
        yline(theta_true(i), '--', 'Color', colors(i,:), 'LineWidth', 2.0); 
    end
    title('Оценка параметров \theta при модифицированном АА');
    xlabel('t'); ylabel('\theta');
    grid on;
    saveas(gcf, 'images/lyon_theta.png');

    % 4. Ошибка компенсации (f - f_hat)
    figure;
    plot(t_lyon, error_f, 'LineWidth', 2);
    title('Ошибка компенсации воздействия f(t) при модифицированном АА');
    xlabel('t'); ylabel('e_f(t)');
    grid on;
    saveas(gcf, 'images/lyon_ef.png');

    % 5. Управление u(t)
    figure;
    plot(t_lyon, u_res, 'LineWidth', 2);
    title('Управление u(t) при модифицированном АА');
    ylim([-1.25 1.5]);
    xlabel('t');
    ylabel('u(t)');
    grid on;
    saveas(gcf, 'images/lyon_u.png');

end

%% Вспомогательные функции

function omega = generate_regressor(t, omega_f, m)
    % omega = [cos(w1*t); sin(w1*t); cos(w2*t); sin(w2*t)]
    num_freqs = length(omega_f);
    omega = zeros(m, 1);
    for i = 1:num_freqs
        omega(2*i-1) = cos(omega_f(i) * t);
        omega(2*i)   = sin(omega_f(i) * t);
    end
end

function dy = sim_lyon(t, y, A, b, K, d, P, omega_f, lambda, gamma, n, m, q, theta_true)
    % Распаковка вектора состояний
    idx = 1;
    x = y(idx : idx + n - 1);           idx = idx + n;
    theta_hat = y(idx : idx + m - 1);   idx = idx + m;
    
    z_filt = reshape(y(idx : idx + n*m - 1), n, m); idx = idx + n*m;
    z_theta_filt = y(idx : idx + n - 1);            idx = idx + n;
    z_omega = y(idx : idx + q*m - 1);               idx = idx + q*m;
    z_omega_theta = y(idx : idx + q - 1);           idx = idx + q;
    z_e = y(idx : idx + q - 1);
    
    % Текущий регрессор
    omega = generate_regressor(t, omega_f, m);
    
    % Истинное возмущение (для моделирования объекта)
    f = theta_true' * omega;
    
    % Управление
    u = -K * x - theta_hat' * omega;

    % 1. Динамика объекта
    dx = A * x + b * u + d * f;
    
    % Вспомогательные величины для фильтров
    Acl = A - b * K;
    
    % 2. Динамика z_filt (фильтрация компонент регрессора через H(s))
    % H(s) = b' * P * (sI - Acl)^-1 * b
    % Реализация: dz/dt = Acl*z + b*u_in; out = b'*P*z
    dz_filt = zeros(n, m);
    omega_f_vec = zeros(m, 1); % Это H(s)[omega]
    
    for j = 1:m
        dz_filt(:, j) = Acl * z_filt(:, j) + b * omega(j);
        omega_f_vec(j) = b' * P * z_filt(:, j);
    end
    
    % 3. Динамика z_theta_filt (фильтрация скаляра omega'*theta_hat)
    val_omega_theta = omega' * theta_hat;
    dz_theta_filt = Acl * z_theta_filt + b * val_omega_theta;
    theta_f_val = b' * P * z_theta_filt; % Это H(s)[omega'*theta_hat]
    
    % 4. Вход для фильтра ошибки e.
    % В схеме Лайона ошибка e - это выходная ошибка слежения.
    % При стабилизации x -> 0, ошибка определяется через b'*P*x
    bPe = b' * P * x; 
    
    % 5. Динамика фильтров DREM (L(s) = lambda / (s+lambda))
    % Применяем к: omega_f_vec (вектор), theta_f_val (скаляр), bPe (скаляр)
    
    dz_omega = zeros(q * m, 1);
    Xi_f = zeros(q, m); % Матрица регрессора для адаптации
    
    % Фильтруем каждый элемент omega_f_vec каждым фильтром lambda_i
    count = 1;
    for i = 1:q
        for j = 1:m
             % Фильтр: dz = -lambda*z + lambda*input
             dz_omega(count) = -lambda(i) * z_omega(count) + lambda(i) * omega_f_vec(j);
             Xi_f(i, j) = z_omega(count);
             count = count + 1;
        end
    end
    
    % Фильтруем theta_f_val -> Xi_f_theta
    dz_omega_theta = zeros(q, 1);
    Xi_f_theta = zeros(q, 1);
    for i = 1:q
        dz_omega_theta(i) = -lambda(i) * z_omega_theta(i) + lambda(i) * theta_f_val;
        Xi_f_theta(i) = z_omega_theta(i);
    end
    
    % Фильтруем bPe -> Xi_e
    dz_e = zeros(q, 1);
    Xi_e = zeros(q, 1);
    for i = 1:q
        dz_e(i) = -lambda(i) * z_e(i) + lambda(i) * bPe;
        Xi_e(i) = z_e(i);
    end
    
    % 6. Закон адаптации (Схема Лайона)
    % dtheta = gamma * Xi_f^T * (Xi_e + Xi_f_theta - Xi_f * theta_hat)
    error_signal = Xi_e + Xi_f_theta - Xi_f * theta_hat;
    dtheta_hat = gamma * Xi_f' * error_signal;
    
    % Сборка производных
    dy = [dx; dtheta_hat; dz_filt(:); dz_theta_filt; dz_omega; dz_omega_theta; dz_e];
end
