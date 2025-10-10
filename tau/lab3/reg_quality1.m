A = [11 -2 13; 6 -1 6; -6 -1 -8];
B = [2; 0; 0];

lambdas = eig(A);

beta = -3;
r = 2;
Q = eye(3);
R = 1;

syms P_ [3, 3]
K_ = -(inv(R + B'*P_*B) * B'*P_*(A - beta*eye(3)));
eqs = (A + B*K_ - beta*eye(3))' * P_ * (A + B*K_ - beta*eye(3)) - r^2*P_ == -Q;
s = vpasolve(eqs, [P_], Random=true);
P = [s.P_1_1, s.P_1_2, s.P_1_3;
     s.P_2_1, s.P_2_2, s.P_2_3;
     s.P_3_1, s.P_3_2, s.P_3_3];
K1 = -double(inv(R + B'*P*B) * B'*P*(A - beta*eye(3)))
ABK1 = A + B*K1 
e1 = eig(A + B*K1)
double(P)

figure;
th = 0:pi/50:2*pi;
xunit = r * cos(th) + beta;
yunit = r * sin(th);
plot(xunit, yunit, "LineWidth", 2);
hold on
plot(real(e1), imag(e1), "o", "LineWidth", 2)
axis equal
grid on
xlabel("Re")
ylabel("Im")
title("Комплексная плоскость и собственные числа при Q = I, R = 1")
hold off
saveas(gcf, "reg_quality1.png")