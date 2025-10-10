% Параметры двигателя
U_1N = 219.393; % Номинальное фазное напряжение, В
r_1 = 1.6; % Активное сопротивление статора, Ом
r_2 = 1.007; % Активное сопротивление ротора, Ом
x_1sigma = 3.0435; % Индуктивное сопротивление рассеяния статора, Ом
x_2sigma = 3.0435; % Индуктивное сопротивление рассеяния ротора, Ом
omega_1 = 314.16; % Угловая частота сети, рад/с
m_1 = 3; % Число фаз
z_p = 2; % Число пар полюсов
k_r = 1.3; % Коэффициент для электромагнитной характеристики (предположение)
k_x = 0.73; % Коэффициент для электромагнитной характеристики (предположение)
% Параметры для расчёта рабочих характеристик
c_1 = 1; % Коэффициент трансформации (предположение)
x_m = 45.988; % Магнитизирующее сопротивление, Ом (из расчёта в отчёте)

% Диапазон скольжения
s = linspace(0.001, 3, 1000); % Избегаем s=0 для избежания деления на ноль

% Механическая характеристика M(s)
M_s = (m_1 * z_p * U_1N^2 * r_2) ./ (omega_1 * s .* ((r_1 + r_2./s).^2 + (x_1sigma + x_2sigma).^2));

% Ток ротора I_2(s)
I_2 = U_1N ./ sqrt((r_1 + r_2./s).^2 + (x_1sigma + x_2sigma).^2);

% Электромагнитная характеристика M_k(s)
M_k_s = (m_1 * z_p * U_1N^2 * r_2 * k_r) ./ (omega_1 * s .* ((r_1 + r_2./s).^2 + (x_1sigma + x_2sigma*k_x).^2));

% Ток ротора I'_2k(s)
I_2k = U_1N ./ sqrt((r_1 + r_2./s).^2 + (x_1sigma + x_2sigma*k_x).^2);

% Рабочие характеристики
P_2 = m_1 .* (I_2.^2) .* r_2 .* (1 - s) ./ s; % Активная мощность на валу (электромагнитная)
I_1 = I_2 + U_1N./(c_1 * x_m); % Приближённо, согласно отчёту
P_1 = P_2 + m_1 .* (I_1.^2) .* r_1 + m_1 .* (I_2.^2) .* r_2; % Подводимая мощность
eta = P_2 ./ P_1; % КПД
cos_phi = P_1 ./ (3 * U_1N .* I_1); % Коэффициент мощности

% Построение графиков — отдельные окна
% Папка для сохранения результатов
figuresDir = fullfile(pwd, 'figures');
if ~exist(figuresDir, 'dir')
mkdir(figuresDir);
end

% График 1: Механическая характеристика M(s)
figure('Name','Механическая характеристика M(s)');
plot(s, M_s, 'b', 'LineWidth', 2);
grid on;
title('Механическая характеристика M(s)');
xlabel('Скольжение s');
ylabel('Момент M, Нм');
% Сохранение графика и данных
exportgraphics(gcf, fullfile(figuresDir, 'M_s.png'), 'Resolution', 300);
exportgraphics(gcf, fullfile(figuresDir, 'M_s.pdf'), 'ContentType', 'vector');
savefig(fullfile(figuresDir, 'M_s.fig'));
writetable(table(s.', M_s.', 'VariableNames', {'s','M'}), fullfile(figuresDir, 'M_s.csv'));

% График 2: Электромагнитная характеристика M_k(s)
figure('Name','Электромагнитная характеристика M_k(s)');
plot(s, M_k_s, 'r', 'LineWidth', 2);
grid on;
title('Электромагнитная характеристика M_k(s)');
xlabel('Скольжение s');
ylabel('Момент M_k, Нм');
% Сохранение графика и данных
exportgraphics(gcf, fullfile(figuresDir, 'M_k_s.png'), 'Resolution', 300);
exportgraphics(gcf, fullfile(figuresDir, 'M_k_s.pdf'), 'ContentType', 'vector');
savefig(fullfile(figuresDir, 'M_k_s.fig'));
writetable(table(s.', M_k_s.', 'VariableNames', {'s','M_k'}), fullfile(figuresDir, 'M_k_s.csv'));

% График 3: Ток I_2(s)
figure('Name','Ток ротора I_2(s)');
plot(s, I_2, 'g', 'LineWidth', 2);
grid on;
title('Ток ротора I_2(s)');
xlabel('Скольжение s');
ylabel('Ток I_2, А');
% Сохранение графика и данных
exportgraphics(gcf, fullfile(figuresDir, 'I_2.png'), 'Resolution', 300);
exportgraphics(gcf, fullfile(figuresDir, 'I_2.pdf'), 'ContentType', 'vector');
savefig(fullfile(figuresDir, 'I_2.fig'));
writetable(table(s.', I_2.', 'VariableNames', {'s','I_2'}), fullfile(figuresDir, 'I_2.csv'));

% График 4: Ток I''_{2k}(s)
figure('Name','Ток ротора I''_{2k}(s)');
plot(s, I_2k, 'm', 'LineWidth', 2);
grid on;
title('Ток ротора I''_{2k}(s)');
xlabel('Скольжение s');
ylabel('Ток I''_{2k}, А');
% Сохранение графика и данных
exportgraphics(gcf, fullfile(figuresDir, 'I_2k.png'), 'Resolution', 300);
exportgraphics(gcf, fullfile(figuresDir, 'I_2k.pdf'), 'ContentType', 'vector');
savefig(fullfile(figuresDir, 'I_2k.fig'));
writetable(table(s.', I_2k.', 'VariableNames', {'s','I_2k'}), fullfile(figuresDir, 'I_2k.csv'));

% График 5: Активная мощность P2(s)
figure('Name','Активная мощность P_2(s)');
plot(s, P_2, 'c', 'LineWidth', 2);
grid on;
title('Активная мощность P_2(s)');
xlabel('Скольжение s');
ylabel('P_2, Вт');
exportgraphics(gcf, fullfile(figuresDir, 'P_2.png'), 'Resolution', 300);
exportgraphics(gcf, fullfile(figuresDir, 'P_2.pdf'), 'ContentType', 'vector');
savefig(fullfile(figuresDir, 'P_2.fig'));
writetable(table(s.', P_2.', 'VariableNames', {'s','P_2'}), fullfile(figuresDir, 'P_2.csv'));

% График 6: Подводимая мощность P1(s)
figure('Name','Подводимая мощность P_1(s)');
plot(s, P_1, 'k', 'LineWidth', 2);
grid on;
title('Подводимая мощность P_1(s)');
xlabel('Скольжение s');
ylabel('P_1, Вт');
exportgraphics(gcf, fullfile(figuresDir, 'P_1.png'), 'Resolution', 300);
exportgraphics(gcf, fullfile(figuresDir, 'P_1.pdf'), 'ContentType', 'vector');
savefig(fullfile(figuresDir, 'P_1.fig'));
writetable(table(s.', P_1.', 'VariableNames', {'s','P_1'}), fullfile(figuresDir, 'P_1.csv'));

% График 7: КПД eta(s)
figure('Name','КПД \eta(s)');
plot(s, eta, 'b--', 'LineWidth', 2);
grid on;
ylim([0 1]);
title('КПД \eta(s)');
xlabel('Скольжение s');
ylabel('\eta');
exportgraphics(gcf, fullfile(figuresDir, 'eta.png'), 'Resolution', 300);
exportgraphics(gcf, fullfile(figuresDir, 'eta.pdf'), 'ContentType', 'vector');
savefig(fullfile(figuresDir, 'eta.fig'));
writetable(table(s.', eta.', 'VariableNames', {'s','eta'}), fullfile(figuresDir, 'eta.csv'));

% График 8: Коэффициент мощности cos(phi)(s)
figure('Name','Коэффициент мощности cos(\phi)(s)');
plot(s, cos_phi, 'r--', 'LineWidth', 2);
grid on;
ylim([0 1]);
title('Коэффициент мощности cos(\phi)(s)');
xlabel('Скольжение s');
ylabel('cos(\phi)');
exportgraphics(gcf, fullfile(figuresDir, 'cos_phi.png'), 'Resolution', 300);
exportgraphics(gcf, fullfile(figuresDir, 'cos_phi.pdf'), 'ContentType', 'vector');
savefig(fullfile(figuresDir, 'cos_phi.fig'));
writetable(table(s.', cos_phi.', 'VariableNames', {'s','cos_phi'}), fullfile(figuresDir, 'cos_phi.csv'));