% ---- параметры системы ----
a = zad1.a;
b = zad1.b;
w = zad1.w;

Td = 0.1;
Tend = 20;

k = 0:Td:Tend;
N = length(k);

% вход
u = sin(w * k);

% выход системы
y = zeros(1, N);
for i = 1:N-1
    y(i+1) = -a * y(i) + b * u(i);
end

% ---- алгоритмы для gamma ----
gammas = [0.5, 10];

figure;

for gi = 1:length(gammas)

    gamma = gammas(gi);

    % параметры оценки
    theta = zeros(2, N);      % [a_hat; b_hat]
    theta_hat = [0; 0];       % начальная оценка

    for k_i = 2:N
        phi = [y(k_i-1); u(k_i-1)];

        % ошибка
        e0 = y(k_i) - phi.'*theta_hat;

        % градиентный алгоритм
        theta_hat = theta_hat + gamma * phi * e0;

        theta(:, k_i) = theta_hat;
    end

    % ---- графики ----
    subplot(2,2,2*gi-1)
    plot(k, theta(1,:), 'LineWidth', 1.6); hold on;
    yline(a, '--'); grid on;
    title(['\gamma = ', num2str(gamma), '  (a\_hat)']);

    subplot(2,2,2*gi)
    plot(k, theta(2,:), 'LineWidth', 1.6); hold on;
    yline(b, '--'); grid on;
    title(['\gamma = ', num2str(gamma), '  (b\_hat)']);
end
