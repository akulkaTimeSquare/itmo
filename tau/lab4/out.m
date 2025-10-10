A = [11 -2 13; 6 -1 6; -6 -1 -8];
B = [2; 0; 0];
C = [2 -2 1];
Bf = [-6 0 0 -1; 0 0 0 0; 6 0 0 0];

G = [35 56 22 -42; -11 -17 -7 12; -6 -10 -5 10; 11 18 6 -13];
Cz = [-2 3 -1];
D = [1 2 1 -1];
Dz = [3 4 2 -3];

lambdas = eig(A);
U1 = [A - lambdas(1)*eye(3) B];
U2 = [A - lambdas(2)*eye(3) B];
U3 = [A - lambdas(3)*eye(3) B];

lambdas(1);
rankU1 = rank(U1);

lambdas(2);
rankU2 = rank(U2);

lambdas(3);
rankU3 = rank(U3);

Q = eye(3);
R = 1;

[P, K, e] = icare(A, B, Q, R);
K1 = -R\B'*P
doubleP = double(P)
er = eig(A + B*K1)

Ahat = [A Bf; zeros(4, 3) G];
Chat = [C D];

Gl = [-5 1 0 0 0 0 0; 0 -5 1 0 0 0 0; 0 0 -5 0 0 0 0; 0 0 0 -6 1 0 0; 0 0 0 0 -6 1 0; 0 0 0 0 0 -6 1; 0 0 0 0 0 0 -7];
Yl = [1; 1; 1; 1; 1; 1; 1];

disp("rank([Yl Gl*Yl Gl^2*Yl Gl^3*Yl Gl^4*Yl Gl^5*Yl Gl^6*Yl])")
rank([Yl Gl*Yl Gl^2*Yl Gl^3*Yl Gl^4*Yl Gl^5*Yl Gl^6*Yl])

Ql = sylvester(Gl, -Ahat, Yl*Chat)
L = Ql\Yl
el = eig(Ahat + L*Chat)

L1 = L(1:3, :)
L2 = L(4:7, :)

syms P [3, 4]
syms Y [1, 4]

eq1 = P*G - A*P == B*Y + Bf;
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

Ad = [A + B*K1 + L1*C, Bf+B*K2+L1*D; L2*C, G + L2*D];
ed = eig(Ad)
eg = eig(G)

function prettyPrintMatrix(A, precision)
    % prettyPrintMatrix(A, precision)
    % Вывод матрицы A без экспоненты
    % precision - число знаков после запятой (по умолчанию 6)

    if nargin < 2
        precision = 3;
    end
    
    fmt = ['%', sprintf('.%df', precision)];
    for i = 1:size(A,1)
        rowStr = '';
        for j = 1:size(A,2)
            rowStr = [rowStr, sprintf([fmt, '\t'], A(i,j))];
        end
        disp(rowStr)
    end
end

disp("P_sol:")
prettyPrintMatrix(double(P_sol), 3)
disp("Y_sol:")
prettyPrintMatrix(double(Y_sol), 3)
for r = 1:size(ed,1)
    for c = 1:size(ed,2)
        fprintf('%8.3f%+8.3fi  ', real(ed(r,c)), imag(ed(r,c)));
    end
    fprintf('\n');
end
