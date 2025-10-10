A = [20 5 -16 9; 6 1 -4 1; 32 9 -25 14; 8 4 -6 4];
C = [-1 0 1 -1];

function prettyPrintMatrix(A, precision)
    % prettyPrintMatrix(A, precision)
    % Вывод матрицы A без экспоненты
    % precision - число знаков после запятой (по умолчанию 6)

    if nargin < 2
        precision = 3;
    end
    
    fmt = ['%', sprintf('.%df', precision)];
    for i = 1:size(A,1)
        rowStr = '';
        for j = 1:size(A,2)
            rowStr = [rowStr, sprintf([fmt, '\t'], A(i,j))];
        end
        disp(rowStr)
    end
end


G1 = [-6 1 0 0; 0 -6 1 0; 0 0 -6 1; 0 0 0 -6];
Y1 = [1; 1; 1; 1];
V1 = [Y1 G1*Y1 G1^2*Y1 G1^3*Y1];
display(rank(V1))
Q1 = sylvester(G1, -A, Y1*C)
L1 = inv(Q1)*Y1
e1 = eig(A + L1*C);
prettyPrintMatrix(A + L1*C)
display(e1)

G2 = [-6 0 0 0; 0 -60 0 0; 0 0 -600 0; 0 0 0 -6000];
Y2 = [1; 1; 1; 1];
V2 = [Y2 G2*Y2 G2^2*Y2 G2^3*Y2];
display(rank(V2))
Q2 = sylvester(G2, -A, Y2*C)
L2 = inv(Q2)*Y2
prettyPrintMatrix(L2);
e2 = eig(A + L2*C);
prettyPrintMatrix(A + L2*C)
prettyPrintMatrix(e2)
display("Alo")

G3 = [-6 7 0 0; -7 -6 0 0; 0 0 -6 8; 0 0 -8 -6];
Y3 = [1; 1; 1; 1];
V3 = [Y3 G3*Y3 G3^2*Y3 G3^3*Y3];
display(rank(V3))
Q3 = sylvester(G3, -A, Y3*C)
L3 = inv(Q3)*Y3;
e3 = eig(A + L3*C);
disp("L3:")
prettyPrintMatrix(L3);
disp("A + L3*C: ")
prettyPrintMatrix(A + L3*C)
display(e3)
