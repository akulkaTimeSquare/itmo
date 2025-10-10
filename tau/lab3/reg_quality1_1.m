A = [11 -2 13; 6 -1 6; -6 -1 -8];
B = [2; 0; 0];
T = [-1 -1.5 -0.5; 0 -1 0; 1 1 0];

hatA = T^(-1)*A*T;
hatB = T^(-1)*B;

hatAc = hatA(2:3, 2:3);
hatBc = hatB(2:3);


lambdas = eig(A);

beta = -3;
r = 2;
Q = eye(2);
R = 1;

syms P_ [2, 2]
K_ = -(inv(R + hatBc'*P_*hatBc) * hatBc'*P_*(hatAc - beta*eye(2)));
eqs = (hatAc + hatBc*K_ - beta*eye(2))' * P_ * (hatAc + hatBc*K_ - beta*eye(2)) - r^2*P_ == -Q;
s = vpasolve(eqs, [P_], Random=true);
P = [s.P_1_1, s.P_1_2;
     s.P_2_1, s.P_2_2];
K = [0 -(inv(R + hatBc'*P*hatBc) * hatBc'*P*(hatAc - beta*eye(2)))]*T^(-1);
e = eig(A + B*K)

figure;
th = 0:pi/50:2*pi;
xunit = r * cos(th) + beta;
yunit = r * sin(th);
plot(xunit, yunit);
hold on
plot(real(e), imag(e), "o")
axis equal
grid on
xlabel("Re")
ylabel("Im")
hold off
saveas(gcf, "reg_quality1_1.png")