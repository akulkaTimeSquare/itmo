omega = logspace(-2, 2, 1000);

denom = sqrt(1 - omega.^2 + omega.^4);

A11 = sqrt(4 + omega.^2) ./ denom;
A12 = 2*sqrt(1 + omega.^2) ./ denom;
A21 = 3*sqrt(1 + omega.^2) ./ denom;
A22 = 6 ./ denom;

L11 = 20*log10(A11);
L12 = 20*log10(A12);
L21 = 20*log10(A21);
L22 = 20*log10(A22);

figure;
semilogx(omega, L11, 'LineWidth', 2); hold on;
semilogx(omega, L12, 'LineWidth', 2);
semilogx(omega, L21, 'LineWidth', 2);
semilogx(omega, L22, 'LineWidth', 2);
hold off;
xlabel('\omega');
ylabel('L(\omega)');
title('ЛАЧХ элементов матрицы');
legend('L_{11}(\omega)', 'L_{12}(\omega)', 'L_{21}(\omega)', 'L_{22}(\omega)');
grid on;
saveas(gcf, 'images/L.png');