if ~exist('images', 'dir')
    mkdir('images');
end

figure;
plot(out.z_feedback, 'LineWidth', 2); hold on;
plot(out.z, 'LineWidth', 2); hold off;
title('Виртуальные выходы z(t) при u = K_1x и u = K_1x + K_2w_g', 'Interpreter', 'tex');
ylabel('z(t)');
xlabel('t');
xlim([0, 5])
grid on;
legend({'z при u = K_1x', 'z при u = K_1x + K_2w_g'});
saveas(gcf, fullfile('images', 'slezh_sost_z_all.png'));

figure;
plot(out.u_feedback, 'LineWidth', 2); hold on;
plot(out.u, 'LineWidth', 2); hold off;
title('Формируемые управления при u = K_1x и u = K_1x + K_2w_g', 'Interpreter', 'tex');
ylabel('u(t)');
xlabel('t');
grid on;
legend({'u = K_1x', 'u = K_1x + K_2w_g'}, 'Location', 'southeast');
saveas(gcf, fullfile('images', 'slezh_sost_u_all.png'));