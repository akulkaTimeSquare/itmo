if ~exist('images', 'dir')
    mkdir('images');
end

figure;
plot(out.z_feedback, 'LineWidth', 2);
title('Виртуальный выход z(t) при u = K_1x', 'Interpreter', 'tex');
ylabel('z(t)');
xlabel('t');
grid on;
saveas(gcf, fullfile('images', 'slezh_sost_z_feedback.png'));

figure;
plot(out.x_feedback, 'LineWidth', 2);
title('Состояния x(t) при u = K_1x', 'Interpreter', 'tex');
ylabel('x(t)');
xlabel('t');
grid on;
legend({'x_1', 'x_2', 'x_3'}, 'Location', 'southeast');
xlim([0 5]);
saveas(gcf, fullfile('images', 'slezh_sost_x_feedback.png'));

figure;
plot(out.u_feedback, 'LineWidth', 2);
title('Формируемое управление u = K_1x', 'Interpreter', 'tex');
ylabel('u(t)');
xlabel('t');
grid on;
saveas(gcf, fullfile('images', 'slezh_sost_u_feedback.png')); 

