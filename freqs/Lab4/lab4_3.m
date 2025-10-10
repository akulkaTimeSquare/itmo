data = readtable("SBER_200411_250410.csv", 'VariableNamingRule', 'preserve');
dates = datetime(string(data.("<DATE>")), 'InputFormat', "yyMMdd");
prices = data.("<CLOSE>");
%plot(prices);

T = 365;
dt = 1;
s = tf('s');
Hs = 1 / (T * s + 1);
Hd = c2d(Hs, dt, 'tustin');
[A, B, C, D] = ssdata(Hd);
x0 = (eye(size(A)) - A) \ B * prices(1);
sys_d = ss(A, B, C, D, dt);

f = figure;
prices_filtered = lsim(sys_d, prices, [], x0);

plot(dates, prices, 'r', 'LineWidth', 2); hold on;
plot(dates, prices_filtered, 'b', 'LineWidth', 2);
legend('Исходные данные', 'Сглаженные данные', 'FontSize', 21);
xlabel('Дата', 'FontSize', 50);
ylabel('Цена', 'FontSize', 50);
title('Сглаживание биржевых данных при T = ' + string(T), 'FontSize', 100);
grid on;

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [100, 100, 1800, 1000]);
set(gca, 'FontSize', 20, 'GridLineWidth', 3, 'LineWidth', 2);
exportgraphics(f, "images\31.jpg")