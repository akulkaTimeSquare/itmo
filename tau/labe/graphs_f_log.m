omega = logspace(-2, 3, 1000);   % частоты

% Знаменатель и его фаза
D = 1 - omega.^2 - 1i*omega;   
phiD = atan2(imag(D), real(D));

% Числители
N11 = 2 - 1i*omega;
N12 = 2*(1 + 1i*omega);
N21 = 3*(1i*omega - 1);
N22 = -6*ones(size(omega));

% Фазы числителей
phiN11 = atan2(imag(N11), real(N11));
phiN12 = atan2(imag(N12), real(N12));
phiN21 = atan2(imag(N21), real(N21));
phiN22 = atan2(imag(N22), real(N22));

% ФЧХ: разность фаз (в градусах)
phi11 = (phiN11 - phiD)*180/pi;
phi12 = (phiN12 - phiD)*180/pi;
phi21 = (phiN21 - phiD)*180/pi;
phi22 = (phiN22 - phiD)*180/pi;

% Построение
figure;
semilogx(omega, phi11, 'LineWidth', 2); hold on;
semilogx(omega, phi12, 'LineWidth', 2);
semilogx(omega, phi21, 'LineWidth', 2);
semilogx(omega, phi22, 'LineWidth', 2);
hold off;
xlabel('\omega'); ylabel('Фаза (град)');
title('ФЧХ элементов матрицы');
legend('f_{11}(\omega)', 'f_{12}(\omega)', 'f_{21}(\omega)', 'f_{22}(\omega)', 'Location', 'northwest');
grid on;
saveas(gcf, 'images/F_log.png');
