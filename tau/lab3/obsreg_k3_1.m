A = [5 -5 -9 3; -5 5 -3 9; -9 -3 5 5; 3 9 5 5];
B = [2; 6; 6; 2];
C = [1 -1 1 1; 1 3 -1 3];

lambdas = eig(A);
ak3 = 4;
R = 1;
nu = 2;
Aa = A + eye(4)*ak3;

Q = e;
[P, K_, e_] = icare(Aa, sqrt(nu)*B, Q, R);
K3 = -inv(R)*B'*P
e = eig(A + B*K3)

function printMatrix(Min, prec)
    fmt = ['%0.', num2str(prec), 'f  '];   % формат вывода
    for i = 1:size(Min,1)
        fprintf(fmt, Min(i,:));
        fprintf('\n');
    end
end

disp("P:")
printMatrix(P, 4)
disp("A + B*K3:")
printMatrix(A+B*K3, 4)