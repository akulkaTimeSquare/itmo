time = out.hatx1.Time;
hatx1 = out.hatx1.Data;
x1 = out.x1.Data;
hatx3 = out.hatx3.Data;
x3 = out.x3.Data;

figure;
plot(time, x1(:,1), 'LineWidth', 2); hold on;
plot(time, hatx1(:,1), 'linestyle', '--', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_1(t)','x̂_1(t)');
title('Первые компоненты векторов состояния при \sigma_1');
saveas(gcf, 'images/x1_1_observer.png');

figure;
plot(time, x1(:,2), 'LineWidth', 2); hold on;
plot(time, hatx1(:,2), 'linestyle', '--', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_2(t)','x̂_2(t)');
title('Вторые компоненты векторов состояния при \sigma_1');
saveas(gcf, 'images/x1_2_observer.png');

figure;
plot(time, x1(:,3), 'LineWidth', 2); hold on;
plot(time, hatx1(:,3), 'linestyle', '--', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_3(t)','x̂_3(t)');
title('Третьи компоненты векторов состояния при \sigma_1');
saveas(gcf, 'images/x1_3_observer.png');

figure;
plot(time, x1(:,4), 'LineWidth', 2); hold on;
plot(time, hatx1(:,4), 'linestyle', '--', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_4(t)','x̂_4(t)');
title('Четвертые компоненты векторов состояния при \sigma_1');
saveas(gcf, 'images/x1_4_observer.png');

figure;
plot(time, x1(:,1) - hatx1(:,1), 'LineWidth', 2); hold on;
plot(time, x1(:,2) - hatx1(:,2), 'LineWidth', 2);
plot(time, x1(:,3) - hatx1(:,3), 'LineWidth', 2);
plot(time, x1(:,4) - hatx1(:,4), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('e_1(t)','e_2(t)','e_3(t)','e_4(t)', 'location', 'southeast');
title('Ошибки наблюдателя e(t)=x(t)-x̂(t) при \sigma_1');
saveas(gcf, 'images/x1_observer_error.png');


figure;
plot(time, x3(:,1), 'LineWidth', 2); hold on;
plot(time, hatx3(:,1), 'linestyle', '--', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_1(t)','x̂_1(t)');
title('Первые компоненты векторов состояния при \sigma_3');
saveas(gcf, 'images/x3_1_observer.png');

figure;
plot(time, x3(:,2), 'LineWidth', 2); hold on;
plot(time, hatx3(:,2), 'linestyle', '--', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_2(t)','x̂_2(t)');
title('Вторые компоненты векторов состояния при \sigma_3');
saveas(gcf, 'images/x3_2_observer.png');

figure;
plot(time, x3(:,3), 'LineWidth', 2); hold on;
plot(time, hatx3(:,3), 'linestyle', '--', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_3(t)','x̂_3(t)');
title('Третьи компоненты векторов состояния при \sigma_3');
saveas(gcf, 'images/x3_3_observer.png');

figure;
plot(time, x3(:,4), 'LineWidth', 2); hold on;
plot(time, hatx3(:,4), 'linestyle', '--', 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('x_4(t)','x̂_4(t)');
title('Четвертые компоненты векторов состояния при \sigma_3');
saveas(gcf, 'images/x3_4_observer.png');

figure;
plot(time, x3(:,1) - hatx3(:,1), 'LineWidth', 2); hold on;
plot(time, x3(:,2) - hatx3(:,2), 'LineWidth', 2);
plot(time, x3(:,3) - hatx3(:,3), 'LineWidth', 2);
plot(time, x3(:,4) - hatx3(:,4), 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('x(t)');
legend('e_1(t)','e_2(t)','e_3(t)','e_4(t)');
title('Ошибки наблюдателя e(t)=x(t)-x̂(t) при \sigma_3');
saveas(gcf, 'images/x3_observer_error.png');