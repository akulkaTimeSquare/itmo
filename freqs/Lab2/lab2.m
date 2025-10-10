T = 100;
N = 16384;
t = linspace(-T, T, N);
dt = t(2) - t(1);
freq = (-N/2:N/2-1)/(N*dt);
omega = 2*pi*freq;

% Наборы параметров [a, b]
params = [1, 1; 1, 4; 2, 1];
colors = ['b', 'r', 'g']; % Цвета для графиков
widths = [13, 9, 4];

%% График всех f(t)
f = figure();
hold on; grid on;
for i = 1:size(params,1)
    a = params(i,1);
    b = params(i,2);
    %f_t = a * double(abs(t) <= b);
    
    %f_t = zeros(size(t));
    %f_t(abs(t) <= b) = a * (1 - abs(t(abs(t) <= b)) / b);
    
    %f_t = a * sinc(b * t/pi);
    
    %f_t = a * exp(-b * t.^2);

    f_t = a * exp(-b * abs(t));
    plot(t, f_t, colors(i), 'LineWidth', widths(i), 'DisplayName', ...
         ['a = ' num2str(a) ', b = ' num2str(b)]);
end
xlabel('t', 'FontSize', 20); ylabel('f(t)', 'FontSize', 20);
xlim([-3.5 3.5]);
ylim([-0.25 2.25]);
title('Функции двустроненного затухания f(t)', 'FontSize', 25);
legend show;
hold off;

%% График всех |F(ω)|
f = figure();
hold on; grid on;
for i = 1:size(params,1)
    a = params(i,1);
    b = params(i,2);
    %F_w = 2*a*b/(sqrt(2*pi))*sinc(omega*b/pi);
    
    %F_w = (a * b / sqrt(2*pi)) * (sinc(omega * b / (2*pi))).^2;
    
    %F_w = (a / b) * sqrt(pi / 2) * double(abs(omega) <= b);

    %F_w = a * sqrt(1 / (2*b)) * exp(-omega.^2 / (4*b));
    
    F_w = (2*a*b) ./ (sqrt(2*pi)*(b^2 + omega.^2));

    plot(omega, F_w, colors(i), 'LineWidth', widths(i), 'DisplayName', ...
         ['a = ' num2str(a) ', b = ' num2str(b)]);
end
xlim([-10 10]);
ylim([-0.15 1.75]);
xlabel('\omega'); ylabel('F(\omega)');
title('Образы Фурье функций двустроненного затухания');
legend show;
hold off;
%plot(nu, hatp, 'r', 'LineWidth', 8);
%xlabel('Частота', "FontSize", 20); ylabel('Амплитуда', 'FontSize', 20);
%title('Фурье-образ П(t) (аналитическое выражение)', 'FontSize', 25);
%ylim([-0.4 1.1]);
%xlim([-5 5]);
%grid on;
%%
a = 2;
b = 1;
%f_t = a * double(abs(t) <= b);
%F_w = 2*a*b/(sqrt(2*pi))*sinc(omega*b/pi);

%f_t = zeros(size(t));
%f_t(abs(t) <= b) = a * (1 - abs(t(abs(t) <= b)) / b);
%F_w = (a * b / sqrt(2*pi)) * (sinc(omega * b / (2*pi))).^2;

%f_t = a * sinc(b * t/pi);
%F_w = (a / b) * sqrt(pi / 2) * double(abs(omega) <= b);

%f_t = a * exp(-b * t.^2);
%F_w = a * sqrt(1 / (2*b)) * exp(-omega.^2 / (4*b));

f_t = a * exp(-b * abs(t));
F_w = (2*a*b) ./ (sqrt(2*pi)*(b^2 + omega.^2));

E_time = trapz(t, abs(f_t).^2)
E_freq = trapz(omega, abs(F_w).^2)
%%
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [100, 100, 1800, 1000]);
set(gca, 'FontSize', 22, 'GridLineWidth', 3, 'LineWidth', 2);
exportgraphics(f, "images\10.jpg")
