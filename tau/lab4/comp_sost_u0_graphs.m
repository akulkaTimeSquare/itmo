if ~exist('images', 'dir')
    mkdir('images');
end

figure;
plot(out.z_u0, 'LineWidth', 2);
title('Виртуальный выход z(t) при u = 0', 'Interpreter', 'tex');
ylabel('z(t)');
xlabel('t');
grid on;
saveas(gcf, fullfile('images', 'comp_sost_z_u0.png'));

figure;
plot(out.x_u0, 'LineWidth', 2);
title('Состояния x(t) при u = 0', 'Interpreter', 'tex');
ylabel('x(t)');
xlabel('t');
grid on;
legend({'x_1', 'x_2', 'x_3'}, 'Location', 'northwest');
saveas(gcf, fullfile('images', 'comp_sost_x_u0.png'));


figure;
plot(out.wf_u0, 'LineWidth', 2);
title('Возмущения w_g(t) при u = 0', 'Interpreter', 'tex');
ylabel('wf(t)');
xlabel('t');
grid on;
legend({'w_{f1}', 'w_{f2}', 'w_{f3}', 'w_{f4}'}, 'Location', 'best');
saveas(gcf, fullfile('images', 'comp_sost_wf_u0.png'));    

