A = [5 -5 -9 3; -5 5 -3 9; -9 -3 5 5; 3 9 5 5];
B = [2; 6; 6; 2];
C = [0 0 0 1; 0 1 0 0];
D = [2; 1];
T = [-1 1 1 -1; 1 1 -1 -1; 1 -1 1 -1; 1 1 1 1];


U = [B, A*B, A*A*B, A*A*A*B];
rankU = rank(U);

V = [C; C*A; C*A*A; C*A*A*A];
rankV = rank(V);

lambdar1 = -13;
lambdar2 = -13;
lambdar3 = -14;
lambdar4 = -14;
Gr = [lambdar1 1 0 0; 0 lambdar2 0 0; 0 0 lambdar3 1; 0 0 0 lambdar4];
Yr = [1 1 1 1];
Pr = sylvester(A, -Gr, B*Yr);
K = -Yr*Pr^(-1);


G = [-5 0; 0 -6];
Y = [-5 0; 1 0];
Ugy = [Y, G*Y];
rankUgy = rank(Ugy);

Q = sylvester(G, -A, Y*C);

CQ = [C; Q];