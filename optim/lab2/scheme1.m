% Временной интервал и сетка
t_span = [0 1];
t_vals = linspace(t_span(1), t_span(2), 201);
init_guess = @(t) [0; 0; 0; 0]; % Начальное приближение [x1; x2; phi1; phi2]
sol_init = bvpinit(t_vals, init_guess);

%% 1. Решение для ОПТИМАЛЬНОГО управления u* (без смещения)
% Формула: u = phi2 / 2
offset_opt = 0; 
sys_opt = @(t, y) ode_system(t, y, offset_opt);
res_bc = @(ya, yb) boundary_conditions(ya, yb);

sol_opt = bvp4c(sys_opt, res_bc, sol_init);
y_opt = deval(sol_opt, t_vals);

x1_opt = y_opt(1, :);
x2_opt = y_opt(2, :);
phi2_opt = y_opt(4, :);
u_opt = phi2_opt / 2 + offset_opt; % Восстанавливаем u по формуле

J_opt = trapz(t_vals, u_opt.^2);
J_cum_opt = cumtrapz(t_vals, u_opt.^2);

fprintf('Оптимальное управление (u*): J = %g\n', J_opt);

%% 2. Решение для управления u1 (смещение -25)
% Мы ищем такую траекторию, где u(t) формируется как phi2/2 - 25,
offset_1 = -25;
sys_1 = @(t, y) ode_system(t, y, offset_1);

sol_1 = bvp4c(sys_1, res_bc, sol_init);
y_1 = deval(sol_1, t_vals);

x1_1 = y_1(1, :);
x2_1 = y_1(2, :);
phi2_1 = y_1(4, :);
u_1 = phi2_1 / 2 + offset_1; 

J_1 = trapz(t_vals, u_1.^2);
J_cum_1 = cumtrapz(t_vals, u_1.^2);

fprintf('Смещенное управление (u* - 25): J = %g\n', J_1);

%% 3. Решение для управления u2 (смещение +50)
% Формула: u = phi2/2 + 50
offset_2 = 50;
sys_2 = @(t, y) ode_system(t, y, offset_2);

sol_2 = bvp4c(sys_2, res_bc, sol_init);
y_2 = deval(sol_2, t_vals);

x1_2 = y_2(1, :);
x2_2 = y_2(2, :);
phi2_2 = y_2(4, :);
u_2 = phi2_2 / 2 + offset_2;

J_2 = trapz(t_vals, u_2.^2);
J_cum_2 = cumtrapz(t_vals, u_2.^2);

fprintf('Смещенное управление (u* + 50): J = %g\n', J_2);

%% x1
figure;
plot(t_vals, x1_opt, 'LineWidth', 2); hold on;
plot(t_vals, x1_1, 'LineWidth', 2);
plot(t_vals, x1_2, 'LineWidth', 2);
hold off;
title('Переменные состояния x_1');
xlabel('t'); ylabel('x_1(t)'); grid on;
legend("u^*", "u_1", "u_2", 'Location', 'northwest');
saveas(gcf, "images/x1m.png");

%% x2
figure;
plot(t_vals, x2_opt, 'LineWidth', 2); hold on;
plot(t_vals, x2_1, 'LineWidth', 2);
plot(t_vals, x2_2, 'LineWidth', 2);
hold off;
title('Переменные состояния x_2');
xlabel('t'); ylabel('x_2(t)'); grid on;
legend("u^*", "u_1", "u_2", 'Location', 'northwest');
saveas(gcf, "images/x2m.png");

%% u
figure;
plot(t_vals, u_opt, 'LineWidth', 2); hold on;
plot(t_vals, u_1, 'LineWidth', 2);
plot(t_vals, u_2, 'LineWidth', 2);
hold off;
title('Формируемые управления');
xlabel('t'); ylabel('u(t)'); grid on;
ylim([-70 105]);
legend("u^*", "u_1", "u_2", 'Location', 'southwest');
saveas(gcf, "images/um.png");

%% J
figure;
plot(t_vals, J_cum_opt, 'LineWidth', 2); hold on;
plot(t_vals, J_cum_1, 'LineWidth', 2);
plot(t_vals, J_cum_2, 'LineWidth', 2);
hold off;
title('Накопленные критерии качества');
xlabel('t'); ylabel('J(t)'); grid on;
legend("J(u^*)", "J(u_1)", "J(u_2)", 'Location', 'northwest');
saveas(gcf, "images/Jm.png");

%% Вспомогательные функции
function dydt = ode_system(t, y, u_offset)
    % y(1)=x1, y(2)=x2, y(3)=phi1, y(4)=phi2
    x1 = y(1);
    x2 = y(2);
    phi1 = y(3);
    phi2 = y(4);
    
    % Управление определяется через сопряженную переменную + смещение
    u = phi2 / 2 + u_offset;
    
    dydt = [
        x2;                         % dx1/dt
        -3*x1 - 6*x2 + u;           % dx2/dt
        3*phi2;                     % dphi1/dt
        -phi1 + 6*phi2              % dphi2/dt
    ];
end

function res = boundary_conditions(ya, yb)
    % Краевые условия:
    % x(0) = [0; 0]
    % x(1) = [10; 0]
    res = [
        ya(1) - 0;
        ya(2) - 0;
        yb(1) - 10;
        yb(2) - 0
    ];
end
