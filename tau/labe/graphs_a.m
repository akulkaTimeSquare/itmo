% Моделирование АЧХ для матрицы
omega = linspace(0, 5, 1000); % диапазон частот

denom = sqrt(1 - omega.^2 + omega.^4);

A11 = sqrt(4 + omega.^2) ./ denom;
A12 = 2*sqrt(1 + omega.^2) ./ denom;
A21 = 3*sqrt(1 + omega.^2) ./ denom;
A22 = 6 ./ denom;

figure;
plot(omega, A11, 'LineWidth', 2); hold on;
plot(omega, A12, 'LineWidth', 2);
plot(omega, A21, 'LineWidth', 2);
plot(omega, A22, 'LineWidth', 2);
xlabel('\omega');
ylabel('A(\omega)');
title('АЧХ элементов матрицы');
legend('|W_{11}(j\omega)|', '|W_{12}(j\omega)|', '|W_{21}(j\omega)|', '|W_{22}(j\omega)|');
grid on;
saveas(gcf, 'images/A.png');