T1 = 0.2;
A1 = -1.3;
omega1 = 0.87;
Am1 = [cos(T1*omega1) sin(T1*omega1); -sin(T1*omega1) cos(T1*omega1)];
H1 = [1 0];
ksi10_1 = 0;
ksi20_1 = A1;
ksi1 = [ksi10_1; ksi20_1];

%% graphs
time = out.ksi.Time;
data = out.ksi.Data;

figure;
stairs(time, data, 'LineWidth', 1.5);
title("Переходный процесс для g(t)=A*sin(omega*t)");
xlabel("t");
ylabel("g(t)");
grid on;
saveas(gcf, 'images/ex_signal_sin.png');