load('ident_lab2_v05.mat');
% исходные параметры
a = zad1.a;
b = zad1.b;
w = zad1.w;

% Параметры
Td = 0.1;
t_end = 10;
t = 0:Td:t_end;

% Входной сигнал
u = sin(w*t);
t1 = 0:0.001:t_end;
u1 = sin(w*t1);

% Создание дискретной передаточной функции
num = [b];
den = [1, a];
Wz = tf(num, den, Td);  % Ts - интервал дискретизации

% Расчет отклика
y = lsim(Wz, u, t);

% Графики
figure;
stairs(t, u, 'LineWidth', 1.5);
hold on;
stairs(t, y, 'LineWidth', 1.5);
hold off;
grid on;
legend('Вход u[k]', 'Выход y[k]');
title('Отклик дискретной системы');
xlabel('t');
ylabel('f(t)');
ylim([-2, 2.15]);
set(gcf, 'PaperPositionMode','auto');
saveas(gcf, "images/task1.png")