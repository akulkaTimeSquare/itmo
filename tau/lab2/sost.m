% out — объект timeseries
t1 = out.x1.Time;      % вектор времени
x1 = out.x1.Data;      % матрица Nx3

figure;
plot(t1, x1(:,1), 'LineWidth', 1.5); hold on;
plot(t1, x1(:,2), 'LineWidth', 1.5);
plot(t1, x1(:,3), 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('x(t)');
xlim([0 6]);
legend('x_1(t)','x_2(t)','x_3(t)');
title('Состояние системы');
yticks(-6:1:7);
ylim([-5.5, 6.5]);
saveas(gcf, 'images/x.png');

figure;
plot(t1, K1_full_standard*x1', 'LineWidth', 1.5);
grid on;
xlabel('t');
xlim([0 6]);
ylabel('u(t)');
title('Управление системой');
saveas(gcf, 'images/u.png');