if ~exist('images', 'dir')
    mkdir('images');
end

figure;
plot(out.z_u0, 'LineWidth', 2);
title('Виртуальный выход z(t) при u = 0', 'Interpreter', 'tex');
ylabel('z(t)');
xlabel('t');
grid on;
saveas(gcf, fullfile('images', 'slezh_sost_z_u0.png'));

figure;
plot(out.x_u0, 'LineWidth', 2);
title('Состояния x(t) при u = 0', 'Interpreter', 'tex');
ylabel('x(t)');
xlabel('t');
grid on;
legend({'x_1', 'x_2', 'x_3'}, 'Location', 'northwest');
saveas(gcf, fullfile('images', 'slezh_sost_x_u0.png'));

figure;
plot(out.wg_u0, 'LineWidth', 2);
title('Возмущения w_g(t) при u = 0', 'Interpreter', 'tex');
ylabel('wg(t)', 'Interpreter', 'tex');
xlabel('t');
grid on;
legend({'w_{g1}', 'w_{g2}', 'w_{g3}', 'w_{g4}'}, 'Location', 'best');
saveas(gcf, fullfile('images', 'slezh_sost_wg_u0.png'));    

