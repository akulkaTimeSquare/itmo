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
plot(time, A1*sin(omega1*time), 'LineWidth', 1.5);
hold on;
stairs(time, data(:, 1), 'LineWidth', 1.5);
hold off;
title("Переходный процесс для g(k) = Asin(kT\omega)", 'Interpreter', 'tex');
xlabel("t");
ylabel("f(t)");
legend("Asin(t\omega)", "g(k)=\xi_1(k)", 'Location', 'northwest');
grid on;
saveas(gcf, 'images/ex_signal_sin.png');