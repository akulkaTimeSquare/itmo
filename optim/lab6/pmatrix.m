A = [3 1; 1 1];
B = [1; 2];
Bf = [3; 4];
Q = [3 0; 0 5];

gamma = 3.25;
P = are(A, B*B'-gamma^(-2)*Bf*Bf', Q);
K = -B'*P;

A_bar = A + B*K;
eig(A_bar)