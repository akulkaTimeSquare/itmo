function printMatrix(Min, prec)
    fmt = ['%0.', num2str(prec), 'f  '];   % формат вывода
    for i = 1:size(Min,1)
        fprintf(fmt, Min(i,:));
        fprintf('\n');
    end
end

A = [20 5 -16 9; 6 1 -4 1; 32 9 -25 14; 8 4 -6 4];
C = [-1 0 1 -1];
x0 = [1; 1; 1; 1];

fa1 = 1/4;
fw1 = 2;

fa2 = 2/4;
fw2 = 4;

fa3 = 3/4;
fw3 = 1;

fa4 = 4/4;
fw4 = 3;

ksia = 1;
ksiw = 5;


lambdas = eig(A);
Q = eye(4);
R = 1;
a = 25;

Q1 = Q;
R1 = R;
[P1, K, e] = icare(A', C', Q1, R1);
L1 = -P1*C'/R1;
e1 = eig(A + L1*C);

disp("P1:")
printMatrix(P1, 4)
disp("L1:")
printMatrix(L1, 4)
disp("A + L1*C:")
printMatrix(A+L1*C, 4)

Q2 = a*Q;
R2 = R;
[P2, K, e] = icare(A', C', Q2, R2);
L2 = -P2*C'/R2;
e2 = eig(A + L2*C);

disp("P2:")
printMatrix(P2, 4)
disp("L2:")
printMatrix(L2, 4)
disp("A + L2*C:")
printMatrix(A+L2*C, 4)

Q3 = Q;
R3 = a*R;
[P3, K, e] = icare(A', C', Q3, R3);
L3 = -P3*C'/R3;
e3 = eig(A + L3*C);

disp("P3:")
printMatrix(P3, 4)
disp("L3:")
printMatrix(L3, 4)
disp("A + L3*C:")
printMatrix(A+L3*C, 4)

Q4 = a*Q;
R4 = a*R;
[P4, K, e] = icare(A', C', Q4, R4);
L4 = -P4*C'/R4;
e4 = eig(A + L4*C);

disp("P4:")
printMatrix(P4, 4)
disp("L4:")
printMatrix(L4, 4)
disp("A + L4*C:")
printMatrix(A+L4*C, 4)

L = L1;