% исходные параметры
a = zad1.a;
b = zad1.b;
w = zad1.w;

Td = 0.1;
Tend = 10;

% дискретное время
k = 0 : Td : Tend;
N = length(k);

% вход
u = sin(w * k);

% выход системы
y = zeros(1, N);

% моделирование y(k+1) = -a*y(k) + b*u(k)
for i = 1:N-1
    y(i+1) = -a * y(i) + b * u(i);
end

% графики
figure; 
plot(k, u, 'LineWidth', 1.5); hold on;
plot(k, y, 'LineWidth', 1.8);
xlabel('t'); grid on;
legend('u(t)', 'y(t)');
title('Modeling of W(z) = b/(z+a)');
