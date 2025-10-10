A = [0 1; -1 1];
B = [1 2; 1 0];
C = [1 -2; 0 3];
D = [-3 0; 0 1];

Bf = [1 2; 1 3];
Df = [1 0; 0 1];
Cz = [1 2; 4 0];
Dz = [4 0; 0 1];

Yg = [0 0 2 0 0 0 0 0;
      0 0 0 -3 0 0 0 0];
Y1 = [0 0 0 0 0 0 7 0;
      0 -3 0 0 0 0 0 0];
Y2 = [0 0 0 0 0 -5 0 0;
      0 0 0 0 0 0 3 0];

x0 = [1; 1];
w0 = [1; 0; 1; 0; 1; 0; 1; 0]

G = [-2 2; -2 -2];
Y = [1 1; 1 1];
Gamma = [0 1 0  0 0  0 0  0;
        -1 0 0  0 0  0 0  0;
        0  0 0  2 0  0 0  0;
        0  0 -2 0 0  0 0  0;
        0  0 0  0 0  3 0  0;
        0  0 0  0 -3 0 0  0;
        0  0 0  0 0  0 0  7;
        0  0 0  0 0  0 -7 0];

P = sylvester(A, -G, B*Y);
disp("P:")
double(P)
Kx = -Y/P;
disp("Kx:")
double(Kx)

lambda = [1i, -1i, 2i, -2i, 3i, -3i, 7i, -7i];
E_lambda1 = [A + B*Kx - eye(2)*lambda(1), B; Cz + Dz*Kx, Dz];
disp("E_lambda1")
double(E_lambda1)
disp("rank(E_lambda1)")
rank(E_lambda1)

E_lambda2 = [A + B*Kx - eye(2)*lambda(2), B; Cz + Dz*Kx, Dz];
disp("E_lambda2")
double(E_lambda2)
disp("rank(E_lambda2)")
rank(E_lambda2)

E_lambda3 = [A + B*Kx - eye(2)*lambda(3), B; Cz + Dz*Kx, Dz];
disp("E_lambda3")
double(E_lambda3)
disp("rank(E_lambda3)")
rank(E_lambda3)

E_lambda4 = [A + B*Kx - eye(2)*lambda(4), B; Cz + Dz*Kx, Dz];
disp("E_lambda4")
double(E_lambda4)
disp("rank(E_lambda4)")
rank(E_lambda4)

E_lambda5 = [A + B*Kx - eye(2)*lambda(5), B; Cz + Dz*Kx, Dz];
disp("E_lambda5")
double(E_lambda5)
disp("rank(E_lambda5)")
rank(E_lambda5)

E_lambda6 = [A + B*Kx - eye(2)*lambda(6), B; Cz + Dz*Kx, Dz];
disp("E_lambda6")
double(E_lambda6)
disp("rank(E_lambda6)")
rank(E_lambda6)

E_lambda7 = [A + B*Kx - eye(2)*lambda(7), B; Cz + Dz*Kx, Dz];
disp("E_lambda7")
double(E_lambda7)
disp("rank(E_lambda7)")
rank(E_lambda7)

E_lambda8 = [A + B*Kx - eye(2)*lambda(8), B; Cz + Dz*Kx, Dz];
disp("E_lambda8")
double(E_lambda8)
disp("rank(E_lambda8)")
rank(E_lambda8)

% Решаем уравнение Франкисона-Дэвисона относительно P и K_omega

% Размерности
n = size(A,1);        % размер состояния
m = size(B,2);        % вход
q = size(Cz,1);       % выход регулирования
r = size(Yg,2);       % число сигналов

% Матрица переменных
P = sym('p', [n r]);   % неизвестная матрица P
Kw = sym('k', [m r]);  % неизвестная матрица Kw
Gw = Gamma;

% Уравнения Франкиса-Дэвисона
eq1 = P*Gw - (A + B*Kx)*P - Bf*Y1 - B*Kw == 0;
eq2 = (Cz + Dz*Kx)*P + Dz*Kw - Yg == 0;

% Разворачиваем
eqs = [eq1(:); eq2(:)];
vars = [P(:); Kw(:)];

% Решение
S = solve(eqs, vars);

% Восстановление решений в матрицы
Psol  = double(reshape(arrayfun(@(x) S.(char(x)), P(:)), n, r));
Kwsol = double(reshape(arrayfun(@(x) S.(char(x)), Kw(:)), m, r));

disp('--- Решение P ---');
disp(Psol);

disp('--- Решение Kw ---');
disp(Kwsol);

ovA = [Gamma zeros(8, 2); Bf*Y1 A];
ovB = [zeros(8, 2); B];
ovC = [Df*Y2-Yg C];
% Проверка наблюдаемости пары (ovC, ovA)
V_ovC_ovA = [ovC; ovC*ovA; ovC*ovA^2; ovC*ovA^3; ovC*ovA^4; ovC*ovA^5; ovC*ovA^6; ovC*ovA^7; ovC*ovA^8]
disp("rank(V_{ovC, ovA})")
rank(V_ovC_ovA)


al = 2;
v = 2;
Q = 0;
R = 1;
nu = 2;
Aa = ovA + eye(size(ovA,1))*al;

[Pl, K_, e_] = icare(Aa', sqrt(nu)*ovC', Q, R);
disp("Pl:")
disp(Pl)
L = -Pl*ovC'*inv(R)
el = eig(ovA + L*ovC)

ovK = [Kwsol Kx]

er = eig(A + B*Kx)

function printMatrix(Min, prec)
    fmt = ['%0.', num2str(prec), 'f  '];   % формат вывода
    for i = 1:size(Min,1)
        fprintf(fmt, Min(i,:));
        fprintf('\n');
    end
end

