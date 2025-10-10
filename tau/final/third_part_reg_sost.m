n = 11;
rng(n , "philox") ;
M = randi([100000 1000000]) / 1000 / sqrt(2) ;
m = randi([1000 10000]) / 1000 * sqrt(3) ;
l = randi([100 1000]) / sqrt(5) / 100;
g = 9.81;

A = [0, 1, 0, 0;
    0, 0, 3*m*g/(4*M + m), 0; 
    0, 0, 0, 1; 
    0, 0, 6*(M + m)*g/l/(4*M + m), 0];

B = [0; 
    4 / (4*M + m);
    0; 
    6/l/(4*M + m)];

C = [1, 0, 0, 0;
    0, 0, 1, 0];

D = [0; 
    6/l/(4*M + m); 
    0; 
    12*(M + m) / m / l^2 / (4*M + m)];

G = [-1  0 0 0; 0 -2 0 0; 0 0 -3 0; 0 0 0 -4];
Y = [1 1 1 1];

rank([Y; Y*G; Y*G^2; Y*G^3])

P = sylvester(A, -G, B*Y);
K = -Y*P^(-1);

e = eig(A + B*K)

function printMatrix(Min, prec)
    fmt = ['%0.', num2str(prec), 'f  '];   % формат вывода
    for i = 1:size(Min,1)
        fprintf(fmt, Min(i,:));
        fprintf('\n');
    end
end

disp("K:");
printMatrix(K, 4);

%% моделирование
function dxdt = nonlinear_pendulum(t, x, u, f, M, m, l, g)
    % x = [x1; x2; x3; x4] = [позиция тележки; скорость тележки; угол маятника; угловая скорость]
    % u - управляющее воздействие на тележку
    % f - внешняя сила на маятник
    
    x1 = x(1); % позиция тележки
    x2 = x(2); % скорость тележки  
    x3 = x(3); % угол маятника
    x4 = x(4); % угловая скорость маятника
    
    % Вычисляем управляющие воздействия (если они функции времени)
    u_val = u;
    f_val = f;
    
    % Нелинейные уравнения движения
    denominator1 = (1/3)*(M + m)*l - (1/4)*m*l*cos(x3)^2;
    numerator1 = (1/3)*l*(u_val - 0.5*m*l*x4^2*sin(x3)) + 0.5*cos(x3)*(f_val + 0.5*m*g*l*sin(x3));
    dx2 = numerator1 / denominator1;
    
    denominator2 = (1/3)*(M + m)*m*l^2 - (1/4)*m^2*l^2*cos(x3)^2;
    numerator2 = 0.5*m*l*cos(x3)*(u_val - 0.5*m*l*x4^2*sin(x3)) + (M + m)*(f_val + 0.5*m*g*l*sin(x3));
    dx4 = numerator2 / denominator2;
    
    dxdt = [x2; dx2; x4; dx4];
end

t_span = [0 10];  % временной интервал
dt = 0.0001;        % шаг интегрирования
t_values = 0:dt:10;

% Начальные условия из CS_project (все состояния ненулевые)
x0a = [0.01;    0;    0;    0];
x0dot_a = [   0; 0.01;    0;    0];
x0phi = [   0;    0; 0.01;    0];
x0dot_phi = [   0;    0;    0; 0.01];

x0a_broken = [20; 0; 0; 0];
x0dot_a_broken = [0; 9; 0; 0];
x0phi_broken = [0; 0; 1.1; 0];
x0dot_phi_broken = [0; 0; 0; 6];

fprintf('Решение a: малые положительные значения [0.01, 0, 0, 0]\n');
[ta, sola] = ode45(@(t,x) nonlinear_pendulum(t, x, K*x, 0, M, m, l, g), t_values, x0a);

fprintf('Решение dot_a: отрицательные значения [0, 0.01, 0, 0]\n');
[tdot_a, soldot_a] = ode45(@(t,x) nonlinear_pendulum(t, x, K*x, 0, M, m, l, g), t_values, x0dot_a);

fprintf('Решение phi: очень малые значения [0, 0, 0.01, 0]\n');
[tphi, solphi] = ode45(@(t,x) nonlinear_pendulum(t, x, K*x, 0, M, m, l, g), t_values, x0phi);

fprintf('Решение dot_phi: отрицательные значения [0, 0, 0, 0.01]\n');
[tdot_phi, soldot_phi] = ode45(@(t,x) nonlinear_pendulum(t, x, K*x, 0, M, m, l, g), t_values, x0dot_phi);

fprintf('Решение a_broken: все состояния ненулевые [20, 0, 0, 0]\n');
[ta_broken, sola_broken] = ode45(@(t,x) nonlinear_pendulum(t, x, K*x, 0, M, m, l, g), 0:dt:2, x0a_broken);

fprintf('Решение dot_a_broken: все состояния ненулевые [0, 9, 0, 0]\n');
[tdot_a_broken, soldot_a_broken] = ode45(@(t,x) nonlinear_pendulum(t, x, K*x, 0, M, m, l, g), 0:dt:2, x0dot_a_broken);

fprintf('Решение phi_broken: все состояния ненулевые [0, 0, 1.1, 0]\n');
[tphi_broken, solphi_broken] = ode45(@(t,x) nonlinear_pendulum(t, x, K*x, 0, M, m, l, g), 0:dt:2, x0phi_broken);

fprintf('Решение dot_phi_broken: все состояния ненулевые [0, 0, 0, 6]\n');
[tdot_phi_broken, soldot_phi_broken] = ode45(@(t,x) nonlinear_pendulum(t, x, K*x, 0, M, m, l, g), 0:dt:3, x0dot_phi_broken);



%% графики
figure;
plot(ta, sola(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(ta, sola(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(ta, sola(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(ta, sola(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('x(t) нелинейной модели при x_0 = [0.01, 0, 0, 0]');
legend('Location', 'northeast');
grid on;
saveas(gcf, 'images/third_part_reg_sost_a.png');

figure;
plot(tdot_a, soldot_a(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(tdot_a, soldot_a(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(tdot_a, soldot_a(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(tdot_a, soldot_a(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('x(t) нелинейной модели при x_0 = [0, 0.01, 0, 0]');
legend('Location', 'northeast');
grid on;
saveas(gcf, 'images/third_part_reg_sost_dot_a.png');

figure;
plot(tphi, solphi(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(tphi, solphi(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(tphi, solphi(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(tphi, solphi(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('x(t) нелинейной модели при x_0 = [0, 0, 0.01, 0]');
legend('Location', 'southeast');
grid on;
saveas(gcf, 'images/third_part_reg_sost_phi.png');

figure;
plot(tdot_phi, soldot_phi(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(tdot_phi, soldot_phi(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(tdot_phi, soldot_phi(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(tdot_phi, soldot_phi(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('x(t) нелинейной модели при x_0 = [0, 0, 0, 0.01]');
legend('Location', 'northeast');
grid on;
saveas(gcf, 'images/third_part_reg_sost_dot_phi.png');

figure;
plot(ta_broken, sola_broken(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(ta_broken, sola_broken(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(ta_broken, sola_broken(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(ta_broken, sola_broken(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('x(t) нелинейной модели при x_0 = [20, 0, 0, 0]');
legend('Location', 'northwest');
grid on;
saveas(gcf, 'images/third_part_reg_sost_a_broken.png');

figure;
plot(tdot_a_broken, soldot_a_broken(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(tdot_a_broken, soldot_a_broken(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(tdot_a_broken, soldot_a_broken(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(tdot_a_broken, soldot_a_broken(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('x(t) нелинейной модели при x_0 = [0, 9, 0, 0]');
legend('Location', 'northwest');
grid on;
saveas(gcf, 'images/third_part_reg_sost_dot_a_broken.png');

figure;
plot(tphi_broken, solphi_broken(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(tphi_broken, solphi_broken(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(tphi_broken, solphi_broken(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(tphi_broken, solphi_broken(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('x(t) нелинейной модели при x_0 = [0, 0, 1.1, 0]');
legend('Location', 'southwest');
grid on;
saveas(gcf, 'images/third_part_reg_sost_phi_broken.png');

figure;
plot(tdot_phi_broken, soldot_phi_broken(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(tdot_phi_broken, soldot_phi_broken(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(tdot_phi_broken, soldot_phi_broken(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(tdot_phi_broken, soldot_phi_broken(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('x(t) нелинейной модели при x_0 = [0, 0, 0, 6]');
legend('Location', 'northwest');
grid on;
saveas(gcf, 'images/third_part_reg_sost_dot_phi_broken.png');
