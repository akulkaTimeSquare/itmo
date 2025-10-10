A = [11 -2 13; 6 -1 6; -6 -1 -8];
B = [2; 0; 0];

G = [35 56 22 -42; -11 -17 -7 12; -6 -10 -5 10; 11 18 6 -13];
Cz = [-2 3 -1];
Dz = [3 4 2 -3];

lambdas = eig(A);
U1 = [A - lambdas(1)*eye(3) B];
U2 = [A - lambdas(2)*eye(3) B];
U3 = [A - lambdas(3)*eye(3) B];

lambdas(1)
rankU1 = rank(U1)

lambdas(2)
rankU2 = rank(U2)

lambdas(3)
rankU3 = rank(U3)

lambdas = eig(A);
Q = zeros(3);
R = 1;

[P, K, e] = icare(A, B, Q, R);
K1 = -inv(R)*B'*P
e1 = eig(A + B*K1)
double(P)

syms P [3, 4]
syms Y [1, 4]

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
P_sol = double(simplify(P_sol));
Y_sol = double(simplify(Y_sol));

function printMatrix(Min, prec)
    fmt = ['%0.', num2str(prec), 'f  '];   % формат вывода
    for i = 1:size(Min,1)
        fprintf(fmt, Min(i,:));
        fprintf('\n');
    end
end

disp("P_sol:")
printMatrix(P_sol, 3)
disp("Y_sol:")
printMatrix(Y_sol, 3)

% если все входные матрицы числовые и нужно число:
% P_num = double(P_sol);  Y_num = double(Y_sol);

K2 = double(Y_sol - K1*P_sol)