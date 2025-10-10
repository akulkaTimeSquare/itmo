A = [5 -5 -9 3; -5 5 -3 9; -9 -3 5 5; 3 9 5 5];
B = [2; 6; 6; 2];
C = [1 -1 1 1; 1 3 -1 3];

lambdas = eig(A);
al1 = 2;
v = 2;
Q = 0;
R = 1;
nu = 2;
Aa = A + eye(4)*al1;

[P, K, e] = icare(Aa', sqrt(nu)*C', Q, R);
L = -P*C'*inv(R)
e = eig(A + L*C)