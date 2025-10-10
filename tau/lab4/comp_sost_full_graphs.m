if ~exist('images', 'dir')
    mkdir('images');
end

figure;
plot(out.z, 'LineWidth', 2);
title('Виртуальный выход z(t) при u = K_1x + K_2w_f', 'Interpreter', 'tex');
ylabel('z(t)');
xlabel('t');
grid on;
saveas(gcf, fullfile('images', 'comp_sost_z_full.png'));

figure;
plot(out.x, 'LineWidth', 2);
title('Состояния x(t) при u = K_1x + K_2w_f', 'Interpreter', 'tex');
ylabel('x(t)');
xlabel('t');
grid on;
legend({'x_1', 'x_2', 'x_3'}, 'Location', 'best');
saveas(gcf, fullfile('images', 'comp_sost_x_full.png'));

figure;
plot(out.u, 'LineWidth', 2);
title('Формируемое управление u = K_1x + K_2w_f', 'Interpreter', 'tex');
ylabel('u(t)');
xlabel('t');
grid on;
saveas(gcf, fullfile('images', 'comp_sost_u_full.png'));