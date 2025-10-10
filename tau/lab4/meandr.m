A = [0 1; 0 0];
B = [0; 1];
C = [1 0];

amplitude = 4;
period = 2*pi;
freq = 2*pi/period;

t = 0:0.01:15;
g = amplitude * square(freq * t);
g_approx = 16/pi * (sin(t) + 1/3*sin(3*t) + 1/5*sin(5*t));

G = [0 1 0 0 0 0; -1 0 0 0 0 0; 0 0 0 3 0 0; 0 0 -3 0 0 0; 0 0 0 0 0 5; 0 0 0 0 -5 0];

Cz = -C;
Dz = 16/pi * [1 0 1/3 0 1/5 0];

figure;
plot(t, g, 'b', t, g_approx, 'r');
grid on;

a1 = 2;
R = 1;
nu = 2;
Aa = A + eye(2)*a1;

Q = zeros(2);
[P, K_, e_] = icare(Aa, sqrt(nu)*B, Q, R);
K1 = -inv(R)*B'*P
disp("P:")
double(P)
e = eig(A + B*K1)

syms P [2, 6]
syms Y [1, 6]

eq1 = P*G - A*P == B*Y;
eq2 = Cz*P + Dz == 0;

% eq1, eq2 заданы как матричные равенства
eqs  = [eq1(:); eq2(:)];          % (добавь eq3(:), если оно тоже есть)
vars = [P(:);  Y(:)];             % неизвестные: все элементы P и Y

% Преобразуем в линейную систему A*x = b и решаем символически
[Aeq, beq] = equationsToMatrix(eqs, vars);
xy = Aeq \ beq;                   % символическое решение

% Возвращаем формы матриц
P_sol = reshape(xy(1:numel(P)), size(P));
Y_sol = reshape(xy(numel(P)+1:end), size(Y));

% по желанию упростить/округлить
P_sol = simplify(P_sol);
Y_sol = simplify(Y_sol);

% если все входные матрицы числовые и нужно число:
% P_num = double(P_sol);  Y_num = double(Y_sol);

K2 = double(Y_sol - K1*P_sol)

function printMatrix(Min, prec)
    fmt = ['%0.', num2str(prec), 'f  '];   % формат вывода
    for i = 1:size(Min,1)
        fprintf(fmt, Min(i,:));
        fprintf('\n');
    end
end

disp("P:")
printMatrix(double(P_sol), 4)
disp("Q:")
printMatrix(double(Y_sol), 4)

x0 = [0; 0];
w0 = [0; 1; 0; 1; 0; 1];


%%
if ~exist('images', 'dir')
    mkdir('images');
end

figure;
plot(out.x, 'LineWidth', 2);
title('Состояния x(t) в задаче слежения за меандром', 'Interpreter', 'tex');
ylabel('x(t)');
xlabel('t');
grid on;
legend({'x_1', 'x_2'}, 'Location', 'northwest');
saveas(gcf, fullfile('images', 'meandr_x.png'));

figure;
plot(out.u, 'LineWidth', 2);
title('Формируемое управление u в задаче слежения за меандром', 'Interpreter', 'tex');
ylabel('u(t)');
xlabel('t');
grid on;
saveas(gcf, fullfile('images', 'meandr_u.png'));

figure;
plot(t, g, 'LineWidth', 2); hold on;
plot(t, g_approx, 'LineWidth', 2);
plot(out.y, 'LineWidth', 2, "Linestyle", "--");
hold off;
title('Графики задающего сигнала, его приближения и выхода', 'Interpreter', 'tex');
ylabel('f(t)');
xlabel('t');
grid on;
legend({'y(t)', 'g(t)', 'g̅(t)'}, 'Location', 'southeast');
saveas(gcf, fullfile('images', 'meandr_y.png')); 