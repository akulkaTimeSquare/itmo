A = [5 -5 -9 3; -5 5 -3 9; -9 -3 5 5; 3 9 5 5];
B = [2; 6; 6; 2];
C = [1 -1 1 1; 1 3 -1 3];

lambdas = eig(A);
al1 = 12;
v = 2;
Q = 0;
R = 1;
nu = 2;
Aa = A + eye(4)*al1;

[P, K_, e] = icare(Aa', sqrt(nu)*C', Q, R);
L1 = -P*C'*inv(R)
e = eig(A + L1*C)
L = L1;


function printMatrix(Min, prec)
    fmt = ['%0.', num2str(prec), 'f  '];   % формат вывода
    for i = 1:size(Min,1)
        fprintf(fmt, Min(i,:));
        fprintf('\n');
    end
end

disp("P:")
printMatrix(P, 4)
disp("A + L1*C:")
printMatrix(A+L1*C, 0)