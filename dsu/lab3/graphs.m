t = out.y.Time;
t_u = 0:T:10;

y = out.y.Data;
u = out.u.Data;

y1 = out.y1.Data;
u1 = out.u1.Data;

y2 = out.y2.Data;
u2 = out.u2.Data;

i = 4;
%% 1
figure;
plot(t, y, 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('y(t)');
title('Выход системы');
ylim([-0.1 1.4]);
xlim([0 10]);
saveas(gcf, 'images/y_' + string(i) + '.png');

figure;
stairs(t_u, u, 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('u(t)');
xlim([0 10]);
ylim([-8 20]);
title('Выход дискретного регулятора');
saveas(gcf, 'images/u_' + string(i) + '.png');

%% 2
figure;
plot(t, y1, 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('y(t)');
title('Выход системы');
ylim([-0.16 0.01]);
xlim([0 10]);
saveas(gcf, 'images/y1_' + string(i) + '.png');

figure;
stairs(t_u, u1, 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('u(t)');
xlim([0 10]);
ylim([-0.1 2.5]);
title('Выход дискретного регулятора');
saveas(gcf, 'images/u1_' + string(i) + '.png');

%% 3
figure;
plot(t, y2, 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('y(t)');
title('Выход системы');
ylim([-0.35 0.35]);
xlim([0 10]);
saveas(gcf, 'images/y2_' + string(i) + '.png');

figure;
stairs(t_u, u2, 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('u(t)');
xlim([0 10]);
ylim([-4.2 5])
title('Выход дискретного регулятора');
saveas(gcf, 'images/u2_' + string(i) + '.png');