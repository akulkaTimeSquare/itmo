B_standard = [2; 0; 0];
A_standard = [11 -2 13; 6 -1 6; -6 -1 -8];

eig(A_standard)

A = [2 3; -3 2];
B = [0; -4];
T = [-1 -1.5 -0.5; 0 -1 0; 1 1 0];

G1 = [-2 1; 0 -2];
Y1 = [1 1];

rank([Y1; Y1*G1]);

P1 = sylvester(A, -G1, B*Y1);
K1 = -Y1*P1^(-1);
K1_full = [2 K1];
K1_full_standard = K1_full*T^(-1);

P1;
K1;
K1_full;
K1_full_standard;


ABK1 = A_standard + B_standard*K1_full_standard;
ABK1;
eig(ABK1);

G2 = [-20 0; 0 -200];
Y2 = [1 1];

rank([Y2; Y2*G2]);

P2 = sylvester(A, -G2, B*Y2);
K2 = -Y2*P2^(-1);
K2_full = [4 K2];
K2_full_standard = K2_full*T^(-1);

P2
K2
K2_full
K2_full_standard


ABK2 = A_standard + B_standard*K2_full_standard;
ABK2;
eig(ABK2);

G3 = [-2 6; -6 -2];
Y3 = [1 1];

rank([Y3; Y3*G3])

P3 = sylvester(A, -G3, B*Y3);
K3 = -Y3*P3^(-1);
K3_full = [6 K3];
K3_full_standard = K3_full*T^(-1);

P3
K3
K3_full
K3_full_standard


ABK3 = A_standard + B_standard*K3_full_standard;
ABK3
eig(ABK3)
