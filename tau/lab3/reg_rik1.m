A = [11 -2 13; 6 -1 6; -6 -1 -8];
B = [2; 0; 0];

lambdas = eig(A);
a1 = 0.5;
R = 1;
nu = 2;
Aa = A + eye(3)*a1;

Q = eye(3);
[P, K_, e_] = icare(Aa, sqrt(nu)*B, Q, R);
K3_1 = -inv(R)*B'*P
e = eig(A + B*K3_1)

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
printMatrix(A+B*K3_1, 4)