K_FB = 1;
K_CO = 7.2;
T = 0.2;

z1 = -1.5;
z2 = -1;
z3 = -0.5;
z4 = 0;
z5 = 0.5;
z6 = 1;
z7 = 1.5;

k1 = (1 - z1)/(K_CO*T);
k2 = (1 - z2)/(K_CO*T);
k3 = (1 - z3)/(K_CO*T);
k4 = (1 - z4)/(K_CO*T);
k5 = (1 - z5)/(K_CO*T);
k6 = (1 - z6)/(K_CO*T);
k7 = (1 - z7)/(K_CO*T);

%% graphs
time = out.x1.Time;
    x1 = out.x1.Data;
x2 = out.x2.Data;
x3 = out.x3.Data;
x4 = out.x4.Data;
x5 = out.x5.Data;
x6 = out.x6.Data;
x7 = out.x7.Data;

function step = step(t)
    step = zeros(size(t));
    step(t >= 1) = 1;
end

figure;
stairs(time, step(time), "LineStyle", '--', 'LineWidth', 1.5);
hold on;
stairs(time, x1, 'LineWidth', 1.5);
hold off;
title("Переходный процесс для z = -1.5");
xlabel("t");
ylabel("f(t)");
grid on;
legend("r(t)", "y(t)", 'Location', 'northwest');
saveas(gcf, 'images/first_case_z1.png');

figure;
stairs(time, step(time), "LineStyle", '--', 'LineWidth', 1.5);
hold on;
stairs(time, x2, 'LineWidth', 1.5);
hold off;
title("Переходный процесс для z = -1");
xlabel("t");
ylabel("f(t)");
grid on;
legend("r(t)", "y(t)", 'Location', 'northwest');
ylim([-0.25 2.25]);
saveas(gcf, 'images/first_case_z2.png');

figure;
stairs(time, step(time), "LineStyle", '--', 'LineWidth', 1.5);
hold on;
stairs(time, x3, 'LineWidth', 1.5);
hold off;
title("Переходный процесс для z = -0.5");
xlabel("t");
ylabel("f(t)");
legend("r(t)", "y(t)");
grid on;
ylim([-0.25 1.6]);
saveas(gcf, 'images/first_case_z3.png');

figure;
stairs(time, step(time), "LineStyle", '--', 'LineWidth', 1.5);
hold on;
stairs(time, x4, 'LineWidth', 1.5);
hold off;
title("Переходный процесс для z = 0");
xlabel("t");
ylabel("f(t)");
legend("r(t)", "y(t)");
grid on;
ylim([-0.25 1.25]);
saveas(gcf, 'images/first_case_z4.png');

figure;
stairs(time, step(time), "LineStyle", '--', 'LineWidth', 1.5);
hold on;
stairs(time, x5, 'LineWidth', 1.5);
hold off;
title("Переходный процесс для z = 0.5");
xlabel("t");
ylabel("f(t)");
grid on;
legend("r(t)", "y(t)");
ylim([-0.25 1.25]);
saveas(gcf, 'images/first_case_z5.png');

figure;
stairs(time, step(time), "LineStyle", '--', 'LineWidth', 1.5);
hold on;
stairs(time, x6, 'LineWidth', 1.5);
hold off;
title("Переходный процесс для z = 1");
xlabel("t");
ylabel("f(t)");
legend("r(t)", "y(t)");
ylim([-0.25 1.25]);
grid on;
saveas(gcf, 'images/first_case_z6.png');

figure;
stairs(time, step(time), "LineStyle", '--', 'LineWidth', 1.5);
hold on;
stairs(time, x7, 'LineWidth', 1.5);
hold off;
title("Переходный процесс для z = 1.5");
xlabel("t");
ylabel("f(t)");
legend("r(t)", "y(t)", 'Location', 'southwest');
grid on;
ylim([-3500 100]);
saveas(gcf, 'images/first_case_z7.png');