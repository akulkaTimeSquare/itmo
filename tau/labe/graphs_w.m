t = linspace(0, 10, 1000);

w11 = exp(t/2).*(-cos(sqrt(3)/2*t) + sqrt(3)*sin(sqrt(3)/2*t));
w12 = exp(t/2).*(2*cos(sqrt(3)/2*t) + 2*sqrt(3)*sin(sqrt(3)/2*t));
w21 = exp(t/2).*(3*cos(sqrt(3)/2*t) - sqrt(3)*sin(sqrt(3)/2*t));
w22 = exp(t/2).*(-4*sqrt(3)*sin(sqrt(3)/2*t));

figure;
plot(t, w11, 'LineWidth', 2); hold on;
plot(t, w12, 'LineWidth', 2);
plot(t, w21, 'LineWidth', 2);
plot(t, w22, 'LineWidth', 2);
hold off;
xlabel('t');
ylabel('w(t)');
title('График весовой характеристики');
legend('w_{11}(t)', 'w_{12}(t)', 'w_{21}(t)', 'w_{22}(t)', "Location", "southwest");
grid on;
saveas(gcf, 'images/W.png');