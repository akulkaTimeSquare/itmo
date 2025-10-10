A = [11 -2 13; 6 -1 6; -6 -1 -8];
B = [2; 0; 0];

lambdas = eig(A);
a1 = 1;
v = 2;
R = 1;
Q = 0;

syms P_ [3, 3]
eqs = [A'*P_ + P_*A - 2*P_'*B*inv(R)*B'*P_ + 2*a1*P_ + Q == 0];
s = vpasolve(eqs, [P_], Random=true);
P = [s.P_1_1, s.P_1_2, s.P_1_3;
     s.P_2_1, s.P_2_2, s.P_2_3;
     s.P_3_1, s.P_3_2, s.P_3_3];
K = -R^(-1)*B'*P
eig(A + B*K)