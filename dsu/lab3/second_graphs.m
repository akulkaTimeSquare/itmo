t = out.x1.Time;
x1 = out.x1.Data;
u1 = out.u1.Data;
x2 = out.x2.Data;
u2 = out.u2.Data;
x3 = out.x3.Data;
u3 = out.u3.Data;
x4 = out.x4.Data;
u4 = out.u4.Data;
x5 = out.x5.Data;
u5 = out.u5.Data;

% График координаты x1 (Случай 1)
figure;
plot(t, x1(:,1), 'LineWidth', 1.5);
hold on;
plot(t, x1(:,2), 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('x(t)');
title('Случай 1: z_1 = 0.1, z_2 = 0.7');
legend('x_1', 'x_2');
saveas(gcf, 'images/first_case.png');