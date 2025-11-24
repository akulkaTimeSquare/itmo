load('ident_lab2_v05.mat');
% исходные параметры
a1 = zad2.a1;
a2 = zad2.a2;
b = zad2.b;
w = zad2.w;

% Параметры
Td = 0.1;
t_end = 15;
t = 0:Td:t_end;
N = t_end/Td+1;

% Входной сигнал
u = sin(w*t);

% Создание дискретной передаточной функции
num = [b];
den = [1, a1, a2];
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
legend('Вход u[k]', 'Выход y[k]', 'Location', 'southeast');
title('Отклик дискретной системы');
xlabel('t');
ylabel('f(t)');
set(gcf, 'PaperPositionMode','auto');
saveas(gcf, "images/task2.png")