% out — объект timeseries
t2 = out.x2.Time;      % вектор времени
x2 = out.x2.Data;      % матрица Nx3

figure;
plot(t2, x2(:,1), 'LineWidth', 1.5); hold on;
plot(t2, x2(:,2), 'LineWidth', 1.5);
plot(t2, x2(:,3), 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_1(t)','x_2(t)','x_3(t)');
title('Состояние системы');
ylim([-4, 2]);
xlim([0 3]);
saveas(gcf, 'images/x1.png');

figure;
plot(t2, K2_full_standard*x2', 'LineWidth', 1.5);
grid on;
xlabel('t');
xlim([0 3]);
ylabel('u(t)');
title('Управление системой');
saveas(gcf, 'images/u1.png');