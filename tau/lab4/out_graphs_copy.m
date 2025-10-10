%% График 1: Траектория управления u(t)
figure('Position', [100, 100, 900, 600]);
plot(out.u, 'LineWidth', 2);
title('Формируемое управление u(t) при z = Cx + Dw = y', 'Interpreter', 'tex', 'FontSize', 14);
ylabel('u(t)', 'FontSize', 12);
xlabel('t', 'FontSize', 12);
grid on;
saveas(gcf, fullfile('images', 'out_u_trajectory.png'));

%% График 2: Ошибка e(t) = x(t) - xhat(t)
figure('Position', [100, 100, 900, 600]);
e1 = out.x - out.xhat;
e2 = out.w - out.what;
plot(e1, 'LineWidth', 2);
hold on;
plot(e2, 'LineWidth', 2, 'LineStyle', '-.');
hold off;
title('Ошибка оценки при z = Cx + Dw = y', 'Interpreter', 'tex', 'FontSize', 14);
ylabel('e(t)', 'FontSize', 12);
xlabel('t', 'FontSize', 12);
grid on;
% Построим корректные подписи легенды для всех кривых e1 и e2
if isprop(e1, 'Data')
    n1 = size(e1.Data, 2);
else
    n1 = size(e1, 2);
end
if isprop(e2, 'Data')
    n2 = size(e2.Data, 2);
else
    n2 = size(e2, 2);
end

legend_labels  = arrayfun(@(i) sprintf('e_{x%d}(t)', i), 1:n1, 'UniformOutput', false);
legend_labels2 = arrayfun(@(i) sprintf('e_{w%d}(t)', i), 1:n2, 'UniformOutput', false);
legend([legend_labels, legend_labels2], 'Location', 'northeast', 'FontSize', 10);
saveas(gcf, fullfile('images', 'out_error_trajectory.png'));

%% График 6: Фазовые траектории y(t) и z(t) на одной координатной плоскости
% Строим фазовые траектории y(t) и z(t) на одной плоскости
figure('Position', [100, 100, 900, 600]);
plot(out.y, 'LineWidth', 2, 'DisplayName', 'y(t)');
hold on;
plot(out.z, 'LineWidth', 2, 'DisplayName', 'z(t)', 'LineStyle', '--');
hold off;

title('Выходы y(t) и z(t) при z = Cx + Dw = y', 'Interpreter', 'tex', 'FontSize', 14);
ylabel('f(t)', 'FontSize', 12);
xlabel('t', 'FontSize', 12);
grid on;
legend('Location', 'northeast', 'FontSize', 12);
saveas(gcf, fullfile('images', 'out_phase_trajectories.png'));

%% График 7: Сравнение состояний x(t) и их оценок xhat(t)
figure('Position', [100, 100, 900, 600]);
t = out.x.Time;
num_states = size(out.x.Data, 2);
hold on;
for i = 1:num_states
    plot(t, out.x.Data(:,i), 'LineWidth', 2);
end
for i = 1:num_states
    plot(t, out.xhat.Data(:,i), 'LineWidth', 2, 'LineStyle', '--');
end
hold off;

legend_x   = arrayfun(@(i) sprintf('x_%d(t)', i), 1:num_states, 'UniformOutput', false);
legend_xhat = arrayfun(@(i) sprintf('x̂_%d(t)', i), 1:num_states, 'UniformOutput', false);
legend([legend_x, legend_xhat], 'Location', 'northeast', 'FontSize', 10);

ylabel('x_i(t)', 'FontSize', 12);
xlabel('t', 'FontSize', 12);
grid on;

title('Сравнение x(t) и x̂(t) при z = Cx + Dw = y', 'Interpreter', 'tex', 'FontSize', 14);
saveas(gcf, fullfile('images', 'out_states_comparison.png'));

%% График 8: Сравнение возмущений w(t) и их оценок what(t)
figure('Position', [100, 100, 900, 600]);
t = out.w.Time;
num_disturbances = size(out.w.Data, 2);
hold on;
for i = 1:num_disturbances
    plot(t, out.w.Data(:,i), 'LineWidth', 2);
end
for i = 1:num_disturbances
    plot(t, out.what.Data(:,i), 'LineWidth', 2, 'LineStyle', '--');
end
hold off;

legend_w   = arrayfun(@(i) sprintf('w_%d(t)', i), 1:num_disturbances, 'UniformOutput', false);
legend_what = arrayfun(@(i) sprintf('ŵ_%d(t)', i), 1:num_disturbances, 'UniformOutput', false);
legend([legend_w, legend_what], 'Location', 'northeast', 'FontSize', 10);

ylabel('w_i(t)', 'FontSize', 12);
xlabel('t', 'FontSize', 12);
grid on;

title('Сравнение w(t) и ŵ(t) при z = Cx + Dw = y', 'Interpreter', 'tex', 'FontSize', 14);
saveas(gcf, fullfile('images', 'out_disturbances_comparison.png'));

