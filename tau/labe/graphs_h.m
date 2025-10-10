t = linspace(0, 10, 1000);

h11 = 2 - 2*exp(t/2).*cos(t*sqrt(3)/2);
h12 = 2 + 2*exp(t/2).*(-cos(t*sqrt(3)/2) + sqrt(3)*sin(t*sqrt(3)/2));
h21 = -3 + exp(t/2).*(3*cos(t*sqrt(3)/2) + sqrt(3)*sin(t*sqrt(3)/2));
h22 = -6 + exp(t/2).*(6*cos(t*sqrt(3)/2 - 2*sqrt(3)*sin(t*sqrt(3)/2)));

figure;
plot(t, h11, 'LineWidth', 2); hold on;
plot(t, h12, 'LineWidth', 2);
plot(t, h21, 'LineWidth', 2);
plot(t, h22, 'LineWidth', 2);
hold off;
xlabel('t');
ylabel('h(t)');
title('График переходной характеристики');
legend('h_{11}(t)', 'h_{12}(t)', 'h_{21}(t)', 'h_{22}(t)', "Location", "northwest");
grid on;
saveas(gcf, 'images/H.png');