n = 11;
rng(n , "philox") ;
M = randi([100000 1000000]) / 1000 / sqrt(2) ;
m = randi([1000 10000]) / 1000 * sqrt(3) ;
l = randi([100 1000]) / sqrt(5) / 100;
g = 9.81;

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

Cz = [0, 0, 1, 0];

% Signal f(t) definition using omega formulation
% f(t) = 0.8*cos(0.5*t + 0.1) + 0.5*cos(1.2*t + 1) + 0.3*cos(2.0*t - 0.5) + 
%        0.2*cos(3.5*t + 0.7) + 0.15*cos(5.0*t - 1.2)

% Frequencies for each sinusoidal component
omega_freqs = [0.5; 1.2; 2.0; 3.5; 5.0];
% Amplitudes for each component
amplitudes = [0.8; 0.5; 0.3; 0.2; 0.15];
% Phase shifts for each component
phases = [0.1; 1; -0.5; 0.7; -1.2];

t = 0:0.0001:10;
f_real = 0.8*cos(0.5*t + 0.1) + 0.5*cos(1.2*t + 1) + 0.3*cos(2.0*t - 0.5) + 0.2*cos(3.5*t + 0.7) + 0.15*cos(5.0*t - 1.2);

G = [   0	0.5	   0	  0	 0	0	   0	  0	 0	0;
-0.5	  0	   0	  0	 0	0	   0	  0	 0	0;
   0	  0	   0	1.2	 0	0	   0	  0	 0	0;
   0	  0	-1.2	  0	 0	0	   0	  0	 0	0;
   0	  0	   0	  0	 0	2	   0	  0	 0	0;
   0	  0	   0	  0	-2	0	   0	  0	 0	0;
   0	  0	   0	  0	 0	0	   0	3.5	 0	0;
   0	  0	   0	  0	 0	0	-3.5	  0	 0	0;
   0	  0	   0	  0	 0	0	   0	  0	 0	5;
   0	  0	   0	  0	 0	0	   0	  0	-5	0];

w0 = [1; 0; 1; 0; 1; 0; 1; 0; 1; 0];
Yf = [0.8*cos(0.1), 0.8*sin(0.1), 0.5*cos(1), 0.5*sin(1), 0.3*cos(-0.5), 0.3*sin(-0.5), 0.2*cos(0.7), 0.2*sin(0.7), 0.15*cos(-1.2), 0.15*sin(-1.2)];

Gr = [-1  0 0 0; 0 -2 0 0; 0 0 -3 0; 0 0 0 -4];
Yr = [1 1 1 1];

P = sylvester(A, -Gr, B*Yr);
K1 = -Yr*P^(-1);

e = eig(A + B*K1)

function printMatrix(Min, prec)
    fmt = ['%0.', num2str(prec), 'f  '];   % формат вывода
    for i = 1:size(Min,1)
        fprintf(fmt, Min(i,:));
        fprintf('\n');
    end
end

disp("K1:");
printMatrix(K1, 4);


syms P [4, 10]
syms Y [1, 10]

eq1 = P*G - A*P == B*Y + D*Yf;
eq2 = Cz*P == 0;

% eq1, eq2 заданы как матричные равенства
eqs  = [eq1(:); eq2(:)];          % (добавь eq3(:), если оно тоже есть)
vars = [P(:);  Y(:)];             % неизвестные: все элементы P и Y

% Преобразуем в линейную систему A*x = b и решаем символически
[Aeq, beq] = equationsToMatrix(eqs, vars);
xy = Aeq \ beq;                   % символическое решение

% Возвращаем формы матриц
P_sol = reshape(xy(1:numel(P)), size(P));
Y_sol = reshape(xy(numel(P)+1:end), size(Y));

K2 = double(Y_sol - K1*P_sol);

disp("K2:");
printMatrix(K2, 0);

x0 = [0.025; -0.06; 0.1; 0.015];
function dxdt = nonlinear_pendulum(t, x, M, m, l, g, K1, K2, Yf, G)
    % x = [x1; x2; x3; x4] = [позиция тележки; скорость тележки; угол маятника; угловая скорость]
    % u - управляющее воздействие на тележку
    % f - внешняя сила на маятник
    
    x1 = x(1); % позиция тележки
    x2 = x(2); % скорость тележки  
    x3 = x(3); % угол маятника
    x4 = x(4); % угловая скорость маятника
    w = x(5:14);
    f_val = Yf*w;
    
    % Вычисляем управляющие воздействия (если они функции времени)
    u_val = K1*[x1; x2; x3; x4] + K2*w;
    
    % Нелинейные уравнения движения
    denominator1 = (1/3)*(M + m)*l - (1/4)*m*l*cos(x3)^2;
    numerator1 = (1/3)*l*(u_val - 0.5*m*l*x4^2*sin(x3)) + 0.5*cos(x3)*(f_val + 0.5*m*g*l*sin(x3));
    dx2 = numerator1 / denominator1;
    
    denominator2 = (1/3)*(M + m)*m*l^2 - (1/4)*m^2*l^2*cos(x3)^2;
    numerator2 = 0.5*m*l*cos(x3)*(u_val - 0.5*m*l*x4^2*sin(x3)) + (M + m)*(f_val + 0.5*m*g*l*sin(x3));
    dx4 = numerator2 / denominator2;

    dw = G*w;
    
    dxdt = [x2; dx2; x4; dx4; dw];
end

function dxdt = linear_pendulum(t, x, A, B, D, K1, K2, Yf, G)
    x_state = x(1:4);
    w = x(5:14);
    
    % Управление
    u = K1*x_state + K2*w;
    
    % Сигнал f
    f = Yf*w;
    
    % Линейная система
    dx_state = A*x_state + B*u + D*f;
    
    % Динамика генератора сигнала
    dw = G*w;
    
    dxdt = [dx_state; dw];
end

t_values = 0:0.0001:10;
[t_nonlin, sol_nonlin] = ode45(@(t,x) nonlinear_pendulum(t, x, M, m, l, g, K1, K2, Yf, G), t_values, [x0; w0]);
u_nonlin = K1*sol_nonlin(:,1:4)' + K2*sol_nonlin(:,5:14)';
x_nonlin = sol_nonlin(:,1:4)';
z_nonlin = Cz*x_nonlin;
w = sol_nonlin(:,5:14)';
f = Yf*w;

[t_lin, sol_lin] = ode45(@(t,x) linear_pendulum(t, x, A, B, D, K1, K2, Yf, G), t_values, [x0; w0]);
u_lin = K1*sol_lin(:,1:4)' + K2*sol_lin(:,5:14)';
x_lin = sol_lin(:,1:4)';
z_lin = Cz*x_lin;
w = sol_lin(:,5:14)';
f = Yf*w;


figure;
plot(t_nonlin, z_nonlin, 'LineWidth', 2);
hold on;
plot(t_lin, z_lin, 'LineWidth', 2, 'LineStyle', '--');
title('Виртуальный выход z(t) при компенсации', 'Interpreter', 'tex');
ylabel('z(t)');
xlabel('t');
grid on;
legend({'Нелинейная система', 'Линейная система'}, 'Location', 'northeast');
saveas(gcf, fullfile('images', 'fifth_part_comp_z.png'));

figure;
plot(t_nonlin, x_nonlin(1, :), 'LineWidth', 2);
hold on;
plot(t_nonlin, x_nonlin(2, :), 'LineWidth', 2);
plot(t_nonlin, x_nonlin(3, :), 'LineWidth', 2);
plot(t_nonlin, x_nonlin(4, :), 'LineWidth', 2);
title('Состояния x(t) нелинейной системы при компенсации', 'Interpreter', 'tex');
ylabel('x(t)');
xlabel('t');
grid on;
legend({'x_1', 'x_2', 'x_3', 'x_4'}, 'Location', 'northwest');
saveas(gcf, fullfile('images', 'fifth_part_comp_x_nonlin.png'));

figure;
plot(t_lin, x_lin(1, :), 'LineWidth', 2);
hold on;
plot(t_lin, x_lin(2, :), 'LineWidth', 2);
plot(t_lin, x_lin(3, :), 'LineWidth', 2);
plot(t_lin, x_lin(4, :), 'LineWidth', 2);
title('Состояния x(t) линейной системы при компенсации', 'Interpreter', 'tex');
ylabel('x(t)');
xlabel('t');
grid on;
legend({'x_1', 'x_2', 'x_3', 'x_4'}, 'Location', 'northwest');
saveas(gcf, fullfile('images', 'fifth_part_comp_x_lin.png'));

figure;
plot(t_nonlin, u_nonlin, 'LineWidth', 2);
hold on;
plot(t_lin, u_lin, 'LineWidth', 2, 'LineStyle', '--');
title('Формируемое управление u = K_1x + K_2w_f', 'Interpreter', 'tex');
ylabel('u(t)');
xlabel('t');
grid on;
legend({'Нелинейная система', 'Линейная система'}, 'Location', 'southeast');
saveas(gcf, fullfile('images', 'fifth_part_comp_u.png'));


figure;
hold on;
plot(t, f_real, 'LineWidth', 2); hold on;
plot(t_nonlin, f, 'LineWidth', 2, 'LineStyle', '--');
title('Сигнал внешнего воздействия f(t)', 'Interpreter', 'tex');
ylabel('f(t)');
xlabel('t');
grid on;
legend({'Истинное воздействие', 'Формируемое генератором'}, 'Location', 'northeast');
saveas(gcf, fullfile('images', 'fifth_part_comp_f.png'));


