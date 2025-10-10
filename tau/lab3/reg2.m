A = [11 -2 13; 6 -1 6; -6 -1 -8];
B = [2; 0; 0];

lambdas = eig(A);
a2 = 1.25;

cvx_begin sdp
variable P(3, 3) symmetric
variable Y(1, 3)
P > 0.0001*eye(3);
P*A' + A*P + 2*a2*P + Y'*B' + B*Y <= 0;
cvx_end

K1_2 = Y*inv(P)

eig(A + B*K1_2)

function printMatrix(Min, prec)
    fmt = ['%0.', num2str(prec), 'f  '];   % формат вывода
    for i = 1:size(Min,1)
        fprintf(fmt, Min(i,:));
        fprintf('\n');
    end
end

disp("P:")
printMatrix(P, 4)
disp("Y:")
printMatrix(Y, 4)
disp("A + B*K1:")
printMatrix(A+B*K1_2, 4)
