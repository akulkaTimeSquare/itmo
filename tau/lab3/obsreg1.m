A = [5 -5 -9 3; -5 5 -3 9; -9 -3 5 5; 3 9 5 5];
B = [2; 6; 6; 2];
C = [1 -1 1 1; 1 3 -1 3];

a1 = 4;
x0 = [1; 1; 1; 1];

cvx_begin sdp
    cvx_solver sdpt3
    cvx_precision high
    variable P(4,4) symmetric
    variable Y(1,4)
    variable mumu nonnegative
    
    minimize(mumu)
    
    % Ограничения
    P > 0.000000000000000000000000000000000001*eye(4);
    P*A' + A*P + 2*a1*P + Y'*B' + B*Y <= 0;
    [P x0;
     x0' 1] > 0;
    [P Y';
     Y mumu] > 0;
cvx_end

K = Y / P   % вместо Y*inv(P)
eig(A + B*K)
