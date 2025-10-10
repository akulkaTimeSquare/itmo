A = [5 -5 -9 3; -5 5 -3 9; -9 -3 5 5; 3 9 5 5];
B = [2; 6; 6; 2];
C = [1 -1 1 1; 1 3 -1 3];

ak3 = 1;
x0 = [1; 1; 1; 1];

cvx_begin sdp
    cvx_precision high
    variable P(4,4) symmetric
    variable Y(1,4)
    variable mumu
    
    minimize(mumu)

    % Ограничения
    P > 0.0001*eye(4);
    P*A' + A*P + 2*ak3*P + Y'*B' + B*Y <= 0;
    [P x0;
     x0' 1] > 0;
    [P Y';
     Y mumu] > 0;
cvx_end

K3 = Y / P   % вместо Y*inv(P)
K = K3;
e = eig(A + B*K3)

function printMatrix(Min, prec)
    fmt = ['%0.', num2str(prec), 'f  '];   % формат вывода
    for i = 1:size(Min,1)
        fprintf(fmt, Min(i,:));
        fprintf('\n');
    end
end

disp("P:")
printMatrix(P, 2)
disp("Y:")
printMatrix(Y, 2)
disp("A + B*K3:")
printMatrix(A+B*K3, 4)
