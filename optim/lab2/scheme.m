%% Решение
syms f1(t) f2(t) x1(t) x2(t)

% Определяем систему ОДУ
ode1 = diff(f1) == 3*f2;
ode2 = diff(f2) == -f1 + 6*f2;
ode3 = diff(x1) == x2;
ode4 = diff(x2) == -3*x1 -6*x2 + f2/2;
odes = [ode1; ode2; ode3; ode4];

% Граничные условия
cond1 = x1(0) == 0;
cond2 = x2(0) == 0;
cond3 = x1(1) == 10;
cond4 = x2(1) == 0;
conds = [cond1; cond2; cond3; cond4];

% Решаем систему
sol = dsolve(odes, conds);

% Извлекаем решения
f1Sol = sol.f1;
f2Sol = sol.f2;
x1Sol = sol.x1;
x2Sol = sol.x2;

% Оптимальное управление
u = f2Sol / 2;

% Упрощаем выражения
f1Sol = simplify(f1Sol);
f2Sol = simplify(f2Sol);
x1Sol = simplify(x1Sol);
x2Sol = simplify(x2Sol);
u = simplify(u);

% Аналитические выражения
disp('f1(t) = '); disp(f1Sol);
disp('f2(t) = '); disp(f2Sol);
disp('x1(t) = '); disp(x1Sol);
disp('x2(t) = '); disp(x2Sol);
disp('u(t) = '); disp(u);

% Вычисляем J
J = int(u^2, t, 0, 1);
J = simplify(J);
J_val = double(J); % численное значение

fprintf('Критерий качества J = ∫₀¹ u(t)² dt = %g\n', J_val);

% Для построения графиков
t_vals = linspace(0, 1, 200);
x1_num = double(subs(x1Sol, t, t_vals));
x2_num = double(subs(x2Sol, t, t_vals));
u_num = double(subs(u, t, t_vals));

%% Графики состояний
figure;
plot(t_vals, x1_num, 'LineWidth', 2); hold on;
plot(t_vals, x2_num, 'LineWidth', 2);
title('Переменные состояния объекта');
xlabel('t');
ylabel('x(t)');
grid on;
legend("x_1(t)", "x_2(t)");
saveas(gcf, "images/x.png");

%% График управления
figure;
plot(t_vals, u_num, 'LineWidth', 2);
title('Оптимальное управление');
xlabel('t');
ylabel('u(t)');
grid on;
saveas(gcf, "images/u.png");

%% Графики накопления
tau = sym('tau');
u_tau = subs(u, t, tau);

% Выражение накопленного критерия
J_cum_sym = int(u_tau^2, tau, 0, t);
J_cum_sym = simplify(J_cum_sym);

% Численная функция для построения графика
J_cum_num = matlabFunction(J_cum_sym, 'Vars', t);

% Вычисляем значения J(t)
u_sq = u_num.^2;                     % u(t)^2 на сетке
J_t_vals = cumtrapz(t_vals, u_sq);  % численный накопленный интеграл

% Проверим, что J(1) = J_val
J_end = J_t_vals(end);
fprintf('Контроль: J(1) = %g, J_val = %g, разница = %g\n', J_end, J_val, abs(J_end - J_val));

%% График накопленного критерия качества J(t)
figure;
plot(t_vals, J_t_vals, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]); % оранжевый
title('Накопленный критерий качества', 'Interpreter', 'latex');
xlabel('t');
ylabel('J(t)');
grid on;
yticks([0 1000 2000 3000 4000 5000 6000 6829])
saveas(gcf, "images/J_cum.png");