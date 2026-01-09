function printMatrix(Min, prec)
    fmt = ['%0.', num2str(prec), 'f  '];   % формат вывода
    for i = 1:size(Min,1)
        fprintf(fmt, Min(i,:));
        fprintf('\n');
    end
end

A = [0 1; 5 -4];
b = [1; 1];
Q = [2 0; 0 1];
r = 5;
x0 = [1; 0];
t_final = 8;
dt = 10e-6;
t = 0:dt:t_final;

[P, K, e] = icare(A, b, Q, r);
A_closed = A - b*K;

jmin = x0'*P*x0;
jexp = 0;
jexp_values = zeros(1, length(t)-1);
x_current = x0;
for i = 1:length(t)-1
    x_t = x_current;
    u_t = -K * x_t;
    integrand = x_t' * Q * x_t + u_t' * r * u_t;
    jexp = jexp + integrand * dt;
    jexp_values(i) = jexp;    
    x_next = x_current + dt * A_closed * x_current;
    x_current = x_next;
end

disp("P:")
printMatrix(P, 3)
disp("K:")
printMatrix(K, 4)
disp("J:")
disp(jexp)


r1 = 15;
k1 = 3;
Q1 = k1*Q;
[P1, K1, e1] = icare(A, b, Q1, r1);
A_closed1 = A - b*K1;

jmin1 = x0'*P1*x0;
jexp1 = 0;
jexp_values1 = zeros(1, length(t)-1);
x_current = x0;
for i = 1:length(t)-1
    x_t = x_current;
    u_t = -K1 * x_t;
    integrand = x_t' * Q1 * x_t + u_t' * r1 * u_t;
    jexp1 = jexp1 + integrand * dt;
    jexp_values1(i) = jexp1;    
    x_next = x_current + dt * A_closed1 * x_current;
    x_current = x_next;
end

disp("P1:")
printMatrix(P1, 3)
disp("K1:")
printMatrix(K1, 4)
disp("J1:")
disp(jexp1)


r2 = 4;
k2 = 2;
Q2 = k2*Q;
[P2, K2, e2] = icare(A, b, Q2, r2);
A_closed2 = A - b*K2;

jmin2 = x0'*P2*x0;
jexp2 = 0;
jexp_values2 = zeros(1, length(t)-1);
x_current = x0;
for i = 1:length(t)-1
    x_t = x_current;
    u_t = -K2 * x_t;
    integrand = x_t' * Q2 * x_t + u_t' * r2 * u_t;
    jexp2 = jexp2 + integrand * dt;
    jexp_values2(i) = jexp2;
    x_next = x_current + dt * A_closed2 * x_current;
    x_current = x_next;
end

disp("P2:")
printMatrix(P2, 3)
disp("K2:")
printMatrix(K2, 4)
disp("J2:")
disp(jexp2)


r3 = 2;
k3 = 4;
Q3 = k3*Q;
[P3, K3, e3] = icare(A, b, Q3, r3);
A_closed3 = A - b*K3;

jmin3 = x0'*P3*x0;
jexp3 = 0;
jexp_values3 = zeros(1, length(t)-1);
x_current = x0;
for i = 1:length(t)-1
    x_t = x_current;
    u_t = -K3 * x_t;
    integrand = x_t' * Q3 * x_t + u_t' * r3 * u_t;
    jexp3 = jexp3 + integrand * dt;
    jexp_values3(i) = jexp3;
    x_next = x_current + dt * A_closed3 * x_current;
    x_current = x_next;
end
disp("P3:")
printMatrix(P3, 3)
disp("K3:")
printMatrix(K3, 4)
disp("J3:")
disp(jexp3)


%% Графики
colors = lines(4);   % хорошая контрастная палитра

figure;
plot(t, jmin * ones(size(t)), 'Color', colors(1,:), 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'J_{min} при (r_0, k_0)');
hold on;
plot(t(1:end-1), jexp_values, 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', 'J_{exp} при (r_0, k_0)');

plot(t, jmin1 * ones(size(t)), 'Color', colors(2,:), 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'J_{min} при (r_1, k_1)');
plot(t(1:end-1), jexp_values1, 'Color', colors(2,:), 'LineWidth', 2, 'DisplayName', 'J_{exp} при (r_1, k_1)');


plot(t, jmin2 * ones(size(t)), 'Color', colors(3,:), 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'J_{min} при (r_2, k_2)');
plot(t(1:end-1), jexp_values2, 'Color', colors(3,:), 'LineWidth', 2, 'DisplayName', 'J_{exp} при (r_2, k_2)');


plot(t, jmin3 * ones(size(t)), 'Color', colors(4,:), 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'J_{min} при (r_3, k_3)');
plot(t(1:end-1), jexp_values3, 'Color', colors(4,:), 'LineWidth', 2, 'DisplayName', 'J_{exp} при (r_3, k_3)');

hold off;
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('J(t)', 'Interpreter', 'tex');
title('Сравнение значений критериев качества', 'Interpreter', 'tex');
legend('Location', 'southeast', 'Interpreter', 'tex');
xlim([0 t_final]);
% saveas(gcf, "images/j_diff.png")