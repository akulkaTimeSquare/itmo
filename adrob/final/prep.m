function printMatrix(Min, prec)
    fmt = ['%0.', num2str(prec), 'f  '];   % формат вывода
    for i = 1:size(Min,1)
        fprintf(fmt, Min(i,:));
        fprintf('\n');
    end
end

A = [0 1; 9 2];
b = [9; 5];
C = [3 3];
d = [9; 5];

Q = eye(2);
r = 3;
[P, K, e] = icare(A, b, Q, r);

disp("P_lqr:")
printMatrix(P, 4)
disp("K:")
printMatrix(K, 4)

Am = A - b*K;
disp("Am:")
printMatrix(Am, 4)
eig(Am)

Q = eye(2);
P = lyap(Am', Q);
disp("P:")
disp(P);

A0f = [0 1 0 0;
       0 0 1 0;
       0 0 0 1;
       -81 -108 -54 -12];
b0f = [0; 0; 0; 1];
N = b0f/d;

eig(A0f)