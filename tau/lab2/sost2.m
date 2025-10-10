% out — объект timeseries
t3 = out.x3.Time;      % вектор времени
x3 = out.x3.Data;      % матрица Nx3

figure;
plot(t3, x3(:,1), 'LineWidth', 1.5); hold on;
plot(t3, x3(:,2), 'LineWidth', 1.5);
plot(t3, x3(:,3), 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('x(t)');
xlim([0 3]);
legend('x_1(t)','x_2(t)','x_3(t)');
title('Состояние системы');
saveas(gcf, 'images/x2.png');

figure;
plot(t3, K3_full_standard*x3', 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('u(t)');
xlim([0 3]);
title('Управление системой');
saveas(gcf, 'images/u2.png');