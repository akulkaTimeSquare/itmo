n = 11;
rng(n , "philox") ;
M = randi([100000 1000000]) / 1000 / sqrt(2) ;
m = randi([1000 10000]) / 1000 * sqrt(3) ;
l = randi([100 1000]) / sqrt(5) / 100;
g = 9.81;

M
m
l
g

A = [0, 1, 0, 0;
    0, 0, 3*m*g/(4*M + m), 0; 
    0, 0, 0, 1; 
    0, 0, 6*(M + m)*g/l/(4*M + m), 0];

B = [0; 
    4 / (4*M + m);
    0; 
    6/l/(4*M + m)];

C = [1, 0, 0, 0;
    0, 0, 1, 0];

D = [0; 
    6/l/(4*M + m); 
    0; 
    12*(M + m) / m / l^2 / (4*M + m)];

A
B
D
C
[V, D_eig] = eig(A);
eigenvalues = diag(D_eig);
eigenvectors = V;
disp('Собственные значения матрицы A:');
disp(eigenvalues);
disp('Собственные векторы матрицы A (столбцы):');
disp(eigenvectors);

U = horzcat(B, A*B, A^2*B, A^3*B)
rank(U)

V = vertcat(C, C*A, C*A^2, C*A^3)
rank(V)

syms s

Wuy = C * inv(s*eye(size(A)) - A) * B;
Wfy = C * inv(s*eye(size(A)) - A) * D;

disp('Полюса и нули передаточной функции Wuy:');
for i = 1:size(Wuy,1)
    for j = 1:size(Wuy,2)
        Wuy_ij = simplify(Wuy(i,j));
        [num, den] = numden(Wuy_ij);
        poles = double(solve(den == 0, s));
        zeros = double(solve(num == 0, s));
        fprintf('Wuy(%d,%d):\n', i, j);
        disp('  Полюса:');
        disp(poles);
        disp('  Нули:');
        disp(zeros);
    end
end

disp('Полюса и нули передаточной функции Wfy:');
for i = 1:size(Wfy,1)
    for j = 1:size(Wfy,2)
        Wfy_ij = simplify(Wfy(i,j));
        [num, den] = numden(Wfy_ij);
        poles = double(solve(den == 0, s));
        zeros = double(solve(num == 0, s));
        fprintf('Wfy(%d,%d):\n', i, j);
        disp('  Полюса:');
        disp(poles);
        disp('  Нули:');
        disp(zeros);
    end
end
