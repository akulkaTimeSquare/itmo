    load('ident_lab2_v05.mat');
a = zad1.a;
b = zad1.b;
w = zad1.w;

Td = 0.1;
t_end = 20;
t = 0:Td:t_end;
N = numel(t);

u = sin(w*t);

Wz = tf([b], [1 a], Td);
y = lsim(Wz, u, t);

gammas = [1, 3, 10];

theta_all = zeros(2, N, numel(gammas));  % theta(2 params, N steps, gamma-index)
err_norms  = zeros(N, numel(gammas));    % ||theta_hat - theta_true||

for gi = 1:numel(gammas)
    gamma = gammas(gi);

    theta_hat = [0; 0];

    for k_i = 2:N
        phi = [-y(k_i-1); u(k_i-1)];
        e0 = y(k_i) - phi.'*theta_hat;
        theta_hat = theta_hat + gamma * (phi * e0)/(1 + gamma*(phi.'*phi));
        theta_all(:, k_i, gi) = theta_hat;

        err = theta_hat - [a; b];
        err_norms(k_i, gi) = norm(err);
    end
end

figure;
subplot(2,1,1); hold on; grid on;
colors = lines(numel(gammas));
for gi = 1:numel(gammas)
    plot(t, squeeze(theta_all(1, :, gi)), 'LineWidth', 1.5, 'Color', colors(gi,:));
end
yline(a, '--k', 'LineWidth', 1.5);
title('Оценки a');
ylim([-0.5 1.1])
yticks([-0.5, 0, 0.5, a]);
xlabel("t")
ylabel("\theta(t)")
legend(arrayfun(@(x) sprintf('\\gamma = %g', x), gammas, 'UniformOutput', false), 'Location', 'southeast');

subplot(2,1,2); hold on; grid on;
for gi = 1:numel(gammas)
    plot(t, squeeze(theta_all(2, :, gi)), 'LineWidth', 1.5, 'Color', colors(gi,:));
end
yline(b, '--k', 'LineWidth', 1.5);
title('Оценки b');
ylim([-0.1 3])
yticks([0 1 2 b])
xlabel("t")
ylabel("\theta(t)")
legend(arrayfun(@(x) sprintf('\\gamma = %g', x), gammas, 'UniformOutput', false), 'Location', 'southeast');

set(gcf, 'PaperPositionMode','auto');
saveas(gcf, "images/task1_all_gammas.png")

% --- Фигура с нормами ошибок ---
figure; hold on; grid on;
for gi = 1:numel(gammas)
    plot(t, err_norms(:, gi), 'LineWidth', 1.5, 'Color', colors(gi,:));
end
title('Нормы параметрических ошибок')
yline(0, '--k', 'LineWidth', 1.5);
legend(arrayfun(@(x) sprintf('\\gamma = %g', x), gammas, 'UniformOutput', false));
ylim([-0.1 3.1]);
xlabel("t")
ylabel("f(t)")
saveas(gcf, "images/task1_error_norms.png")
