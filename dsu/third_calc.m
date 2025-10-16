T = 0.2;
A = 5;
B = 5.5;
C = 1.5;
Am = [0 1 0; 0 0 1; 1 -3 3];
ksi10 = A;
ksi20 = A + B*T + C*T^2;
ksi30 = A + 2*B*T + 3*C*T^2;
ksi = [ksi10; ksi20; ksi30];

%% graphs
time = out.ksi1.Time;
data = out.ksi1.Data;

figure;
stairs(time, data, 'LineWidth', 1.5);
title("Переходный процесс для g(t)=A+BkT+CkT^2");
xlabel("t");
ylabel("g(t)");
grid on;
saveas(gcf, 'images/ex_signal.png');