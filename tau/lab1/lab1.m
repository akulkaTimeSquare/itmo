
A = [7, -6, 9; 6, -5, 6; -6, 3, -8];
B = [-2; -1; 2];
x1 = [-5; -3; 3];
t1 = 3;

x1 = [-5; -3;  3];
P_t1 = [0.3688436, 0.20346076, -0.37345989;
        0.20346076, 0.13461396, -0.17461489;
        -0.37345989, -0.17461489, 0.46461429];

t1 = 3;

P_inv = inv(P);

syms t;
exp_A_transpose = expm(A' * (t1 - t));  % e^{A^T (t1 - t)}
u_t = B' * exp_A_transpose * inv(P_t1) * x1;

disp(simplify(expm(A*t)*B*B'*expm(A'*t)));
%u_t_simplified = simplify(u_t);
%disp('Управление u(t):');
%disp(u_t_simplified);

% Определяем функцию правых частей системы
ode_fun = @(t, x) A * x + B * double(subs(u_t, t));

% Начальные условия (x(0) = 0, так как управление переводит систему из 0 в x1)
x0 = [0; 0; 0];

% Временной интервал
t_span = [0, t1];

% Решаем ОДУ
[t_sol, x_sol] = ode45(ode_fun, t_span, x0);

% Проверяем, что x(t1) ≈ x1
x_final = x_sol(end, :)';
disp('x(t1):');
disp(x_final);
disp('Ожидаемое x1:');
disp(x1);

figure;
plot(t_sol, x_sol(:, 1), 'r', t_sol, x_sol(:, 2), 'g', t_sol, x_sol(:, 3), 'b');
legend('x_1(t)', 'x_2(t)', 'x_3(t)');
xlabel('Время t');
ylabel('Состояние x(t)');
title('Траектория системы под управлением u(t)');
grid on;
