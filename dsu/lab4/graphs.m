t = out.y.Time;
t_u = 0:1:40;

y = out.y.Data;
u = out.u.Data;

y1 = out.y1.Data;
u1 = out.u1.Data;
%% 0
figure;
stairs(t, y, 'LineWidth', 2); hold on;
stairs(t, double(t >= 1), 'LineWidth', 2, 'LineStyle', '--');
hold off;
grid on;
xlabel('t');
ylabel('y(t)');
title('Выход системы при апериодическом регуляторе');
xlim([0 10]);
ylim([-0.05 1.05]);
legend("Выход системы", "Задающее воздействие", 'Location','southeast');
saveas(gcf, 'images/y.png');

figure;
stairs(t_u, u, 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('u(t)');
xlim([0 10]);
ylim([-0.5 7.5]);
title('Управляющее воздействие при апериодическом регуляторе');
saveas(gcf, 'images/u.png');

%% 1
figure;
stairs(t, y1, 'LineWidth', 2); hold on;
stairs(t, double(t >= 1), 'LineWidth', 2, 'LineStyle', '--');
hold off;
grid on;
xlabel('t');
ylabel('y(t)');
title('Выход системы при регуляторе Далина');
legend("Выход системы", "Задающее воздействие", 'Location','southeast');
xlim([0 40]);
ylim([-0.05 1.05]);
saveas(gcf, 'images/y1.png');

figure;
stairs(t_u, u1, 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('u(t)');
xlim([0 10]);
ylim([-0.05 1.05]);
title('Управляющее воздействие при регуляторе Далина');
saveas(gcf, 'images/u1.png');