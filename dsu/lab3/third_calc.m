T = 0.25;
A = 5;
B = 5.5;
C = 1.5;
Am = [0 1 0; 0 0 1; 1 -3 3];
ksi10 = A;
ksi20 = A + B*T + C*T^2;
ksi30 = A + 2*B*T + 4*C*T^2;
ksi = [ksi10; ksi20; ksi30];

%% graphs
time = out.ksi1.Time;
data = out.ksi1.Data;

figure;
plot(time, A + B*time + C*time.^2, 'LineWidth', 1.5);
hold on;
stairs(time, data, 'LineWidth', 1.5);
hold off;
title("Переходный процесс для g[k]=A+BkT+C(kT)^2");
xlabel("t");
ylabel("f(t)");
grid on;
xlim([0 10]);
legend("A+Bt+Ct^2", "g(k)=\xi_1(k)", 'Location', 'northwest');
saveas(gcf, 'images/ex_signal.png');