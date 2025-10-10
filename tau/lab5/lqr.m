function printMatrix(Min, prec)
    fmt = ['%0.', num2str(prec), 'f  '];   % формат вывода
    for i = 1:size(Min,1)
        fprintf(fmt, Min(i,:));
        fprintf('\n');
    end
end

A = [11 -2 13; 6 -1 6; -6 -1 -8];
B = [2; 0; 0];
x0 = [1; 1; 1];

lambdas = eig(A);
Q = eye(3);
R = 1;
a = 15;


Q1 = Q;
R1 = R;
[P1, K, e] = icare(A, B, Q1, R1);
K1 = -inv(R1)*B'*P1;
e1 = eig(A + B*K1);

disp("P1:")
printMatrix(P1, 4)
disp("K1:")
printMatrix(K1, 4)
disp("A + BK1:")
printMatrix(A+B*K1, 4)


Q2 = a*Q;
R2 = R;
[P2, K, e] = icare(A, B, Q2, R2);
K2 = -inv(R2)*B'*P2;
e2 = eig(A + B*K2);

disp("P2:")
printMatrix(P2, 4)
disp("K2:")
printMatrix(K2, 4)
disp("A + BK2:")
printMatrix(A+B*K2, 4)

Q3 = Q;
R3 = a*R;
[P3, K, e] = icare(A, B, Q3, R3);
K3 = -inv(R3)*B'*P3;
e3 = eig(A + B*K3);

disp("P3:")
printMatrix(P3, 4)
disp("K3:")
printMatrix(K3, 4)
disp("A + BK3:")
printMatrix(A+B*K3, 4)

Q4 = a*Q;
R4 = a*R;
[P4, K, e] = icare(A, B, Q4, R4);
K4 = -inv(R4)*B'*P4;
e4 = eig(A + B*K4);

disp("P4:")
printMatrix(P4, 4)
disp("K4:")
printMatrix(K4, 4)
disp("A + BK4:")
printMatrix(A+B*K4, 4)

A_closed1 = A + B*K1;
A_closed2 = A + B*K2;
A_closed3 = A + B*K3;
A_closed4 = A + B*K4;

t_final = 2;
dt = 10e-6;
t = 0:dt:t_final;

jmin1 = x0'*P1*x0;
jexp1 = 0;
jexp_values1 = zeros(1, length(t)-1);
x_current = x0;
for i = 1:length(t)-1
    x_t = x_current;
    u_t = K1 * x_t;
    integrand = x_t' * Q1 * x_t + u_t' * R1 * u_t;
    jexp1 = jexp1 + integrand * dt;
    jexp_values1(i) = jexp1;    
    x_next = x_current + dt * A_closed1 * x_current;
    x_current = x_next;
end

jmin2 = x0'*P2*x0;
jexp2 = 0;
jexp_values2 = zeros(1, length(t)-1);
x_current = x0;
for i = 1:length(t)-1
    x_t = x_current;
    u_t = K2 * x_t;
    integrand = x_t' * Q2 * x_t + u_t' * R2 * u_t;
    jexp2 = jexp2 + integrand * dt;
    jexp_values2(i) = jexp2;    
    x_next = x_current + dt * A_closed2 * x_current;
    x_current = x_next;
end

jmin3 = x0'*P3*x0;
jexp3 = 0;
jexp_values3 = zeros(1, length(t)-1);
x_current = x0;
for i = 1:length(t)-1
    x_t = x_current;
    u_t = K3 * x_t;
    integrand = x_t' * Q3 * x_t + u_t' * R3 * u_t;
    jexp3 = jexp3 + integrand * dt;
    jexp_values3(i) = jexp3;    
    x_next = x_current + dt * A_closed3 * x_current;
    x_current = x_next;
end

jmin4 = x0'*P4*x0;
jexp4 = 0;
jexp_values4 = zeros(1, length(t)-1);
x_current = x0;
for i = 1:length(t)-1
    x_t = x_current;
    u_t = K4 * x_t;
    integrand = x_t' * Q4 * x_t + u_t' * R4 * u_t;
    jexp4 = jexp4 + integrand * dt;
    jexp_values4(i) = jexp4;    
    x_next = x_current + dt * A_closed4 * x_current;
    x_current = x_next;
end


%% Построение графика jexp(t)
figure;
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.05 0.05 0.4 0.5]);
plot(t(1:end-1), jexp_values1, 'r-', 'LineWidth', 2, 'DisplayName', 'J_{exp}(t) при (Q,R)');
hold on;
plot(t, jmin1 * ones(size(t)), 'r--', 'LineWidth', 2, 'DisplayName', 'J_{min} при (Q,R)');
hold on;
plot(t(1:end-1), jexp_values2, 'g-', 'LineWidth', 2, 'DisplayName', 'J_{exp}(t) при (aQ,R)');
hold on;
plot(t, jmin2 * ones(size(t)), 'g--', 'LineWidth', 2, 'DisplayName', 'J_{min} при (aQ,R)');
hold on;
plot(t(1:end-1), jexp_values3, 'b-', 'LineWidth', 2, 'DisplayName', 'J_{exp}(t) при (Q,aR)');
hold on;
plot(t, jmin3 * ones(size(t)), 'b--', 'LineWidth', 2, 'DisplayName', 'J_{min} при (Q,aR)');
hold on;
plot(t(1:end-1), jexp_values4, 'm-', 'LineWidth', 2, 'DisplayName', 'J_{exp}(t) при (aQ,aR)');
hold on;
plot(t, jmin4 * ones(size(t)), 'm--', 'LineWidth', 2, 'DisplayName', 'J_{min} при (aQ,aR)');
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('J(t)', 'Interpreter', 'tex');
title('Функционалы качества J(t)', 'Interpreter', 'tex', 'FontSize', 13);
legend('Location', 'best', 'Interpreter', 'tex');
xlim([0 t_final]);

% Добавление отступа сверху и значений jmin1..jmin4 в метки оси Y
ax = gca;
ymax = max([jexp_values1(end), jexp_values2(end), jexp_values3(end), jexp_values4(end), jmin1, jmin2, jmin3, jmin4]);
ymin = 0;
margin = 0.12 * max(1, ymax - ymin);
ylim([ymin, ymax + margin]);
yt = get(ax, 'YTick');
yt = unique([yt, jmin1, jmin2, jmin3, jmin4]);
set(ax, 'YTick', sort(yt));

% Сохранение графика
%saveas(gcf, 'images/jexp_plot.png');