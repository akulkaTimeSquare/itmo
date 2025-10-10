A = [5 -5 -9 3; -5 5 -3 9; -9 -3 5 5; 3 9 5 5];
B = [2; 6; 6; 2];
C = [1 -1 1 1; 1 3 -1 3];
D = [2; 1];
T = [-1 1 1 -1; 1 1 -1 -1; 1 -1 1 -1; 1 1 1 1];

hatA = inv(T)*A*T;
hatB = inv(T)*B;
hatC = C*T;

hatAc = hatA(1:3, 1:3);
hatBc = hatB(1:3, 1);
hatCc = hatC(1:2, 1:3);

U = [B, A*B, A*A*B, A*A*A*B];
rankU = rank(U);

V = [C; C*A; C*A*A; C*A*A*A];
rankV = rank(V);

lambdas = eig(A);
V1 = [A - lambdas(1) * eye(4); C];
rankV1 = rank(V1);
V2 = [A - lambdas(2) * eye(4); C];
rankV2 = rank(V2);
V3 = [A - lambdas(3) * eye(4); C];
rankV3 = rank(V3);
V4 = [A - lambdas(4) * eye(4); C];
rankV4 = rank(V4);

lambdar1 = -13;
lambdar2 = -13;
lambdar3 = -14;
lambdar4 = -14;

lambdan1 = -11;
lambdan2 = -11;
lambdan3 = -12;

Gr = [lambdar1 1 0 0; 0 lambdar2 0 0; 0 0 lambdar3 1; 0 0 0 lambdar4];
Gn = [lambdan1 1 0; 0 lambdan2 0; 0 0 lambdan3];

Yr = [1 1 1 1];
Vr = [Yr; Yr*Gr; Yr*Gr*Gr; Yr*Gr*Gr*Gr];
rankVr = rank(Vr);

Yn = [1 1; 1 1; 1 1];
Un = [Yn, Gn*Yn, Gn*Gn*Yn];
rankUn = rank(Un);

Pr = sylvester(A, -Gr, B*Yr);
K = -Yr*Pr^(-1);

Qn = sylvester(Gn, -hatAc, Yn*hatCc);
hatLnc = Qn^(-1)*Yn;
hatLn = [hatLnc; 0 0];
L = T*hatLn;

Ar = A + B*K;
lambdasAr = eig(Ar)

An = A + L*C;
eig(An)

function printMatrix(Min, prec)
    fmt = ['%0.', num2str(prec), 'f  '];   % формат вывода
    for i = 1:size(Min,1)
        fprintf(fmt, Min(i,:));
        fprintf('\n');
    end
end

disp("Ar:")
printMatrix(Ar, 11)
disp("Q:")
printMatrix(Qn, 4)
disp("hatLnc")
printMatrix(hatLnc, 3)
disp("L")
printMatrix(L, 4)
disp("A + LC")
printMatrix(A+L*C, 4)
disp("Pr:")
printMatrix(Pr, 4)
disp("K:")
printMatrix(K, 4)