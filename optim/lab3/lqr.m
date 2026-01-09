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

[P, K, e] = icare(A, b, Q, r);

disp("P:")
printMatrix(P, 3)
disp("K:")
printMatrix(K, 4)
A_closed = A - b*K;

t_final = 8;
dt = 10e-6;
t = 0:dt:t_final;

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

%% Графики
figure;
plot(t, jmin * ones(size(t)), 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'J_{min}');
hold on;
plot(t(1:end-1), jexp_values, 'LineWidth', 2, 'DisplayName', 'J_{exp}');
hold off;
grid on;
xlabel('t', 'Interpreter', 'tex');
ylabel('J(t)', 'Interpreter', 'tex');
title('Критерий качества', 'Interpreter', 'tex');
legend('Location', 'southeast', 'Interpreter', 'tex');
yticks([0 1 2 3 4 5 6 7 7.8]);
xlim([0 t_final]);
saveas(gcf, "images/j.png")