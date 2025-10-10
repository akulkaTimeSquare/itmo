time = out.x2.Time;
hatx2 = out.hatx2.Data;
x2 = out.x2.Data;

figure;
plot(time, x2(:,1), 'LineWidth', 2); hold on;
plot(time, hatx2(:,1), 'linestyle', '-.', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_1(t)','x̂_1(t)', 'location', 'southeast');
title('Первые компоненты векторов состояния при \sigma_2');
saveas(gcf, 'images/x2_1_observer.png');

figure;
plot(time, x2(:,2), 'LineWidth', 2); hold on;
plot(time, hatx2(:,2), 'linestyle', '-.', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_2(t)','x̂_2(t)', 'location', 'northeast');
title('Вторые компоненты векторов состояния при \sigma_2');
saveas(gcf, 'images/x2_2_observer.png');

figure;
plot(time, x2(:,3), 'LineWidth', 2); hold on;
plot(time, hatx2(:,3), 'linestyle', '-.', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_3(t)','x̂_3(t)', 'location', 'southeast');
title('Третьи компоненты векторов состояния при \sigma_2');
saveas(gcf, 'images/x2_3_observer.png');

figure;
plot(time, x2(:,4), 'LineWidth', 2); hold on;
plot(time, hatx2(:,4), 'linestyle', '-.', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_4(t)','x̂_4(t)', 'location', 'southeast');
title('Четвертые компоненты векторов состояния при \sigma_2');
saveas(gcf, 'images/x2_4_observer.png');

figure;
plot(time, x2(:,1) - hatx2(:,1), 'LineWidth', 2); hold on;
plot(time, x2(:,2) - hatx2(:,2), 'LineWidth', 2);
plot(time, x2(:,3) - hatx2(:,3), 'LineWidth', 2);
plot(time, x2(:,4) - hatx2(:,4), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('e_1(t)','e_2(t)','e_3(t)','e_4(t)');
title('Ошибки наблюдателя e(t)=x(t)-x̂(t) при \sigma_2');
saveas(gcf, 'images/x2_observer_error.png');