A = [5 -5 -9 3; -5 5 -3 9; -9 -3 5 5; 3 9 5 5];
B = [2 0; 6 0; 6 0; 2 0];
C = [1 -1 1 1; 1 3 -1 3];
D = [0 2; 0 1];

x0 = [1; 1; 1; 1];
xhat0 = [0; 0; 0; 0];

F = [1 0 0 0; 0 2 0 0; 0 0 4 0; 0 0 0 3];
E = [2 0; 0 5];

Qk = eye(4);
Rk = eye(2);

[Pk, Kk, ek] = icare(A, B, Qk, Rk);
K = -Rk\B'*Pk;
e = eig(A + B*K);

Ql = F;
Rl = E;

[Pl, Kl, el] = icare(A', C', Ql, Rl);
L = -Pl*C'/Rl;

function printMatrix(Min, prec)
    fmt = ['%0.', num2str(prec), 'f  '];   % формат вывода
    for i = 1:size(Min,1)
        fprintf(fmt, Min(i,:));
        fprintf('\n');
    end
end

disp("Pk:")
printMatrix(Pk, 4)
disp("K:")
printMatrix(K, 4)
disp("Pl:")
printMatrix(Pk, 4)
disp("L:")
printMatrix(L, 4)