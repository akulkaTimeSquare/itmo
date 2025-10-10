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

x0_1 = [0.01; 0.01; 0.01; 0.01];
x0_2 = [0.1;  0; 0.01;  0];

a1 = 0.5;
a2 = 1.5;

cvx_begin sdp
    variable P1(4, 4) symmetric
    variable Y1(1, 4)
    variable mumu1
    minimize mumu1
    P1 > 0.0001*eye(4);
    P1*A' + A*P1 + 2*a1*P1 + Y1'*B' + B*Y1 <= 0;
    [P1 Y1';
    Y1 mumu1] > 0;
    [P1 x0_1;
    x0_1' 1] > 0;
cvx_end

cvx_begin sdp
    variable P2(4,4) symmetric
    variable Y2(1,4)
    variable mumu2
    minimize mumu2
    P2 > 1e-4*eye(4);
    P2*A' + A*P2 + 2*a2*P2 + Y2'*B' + B*Y2 <= 0;
    [P2 Y2'; Y2 mumu2] > 0;
    [P2 x0_2; x0_2' 1] > 0;
cvx_end

K1 = Y1/P1;
K2 = Y2/P2;
e1 = eig(A + B*K1)
e2 = eig(A + B*K2)

%x0_2_1 = [0.1;  0; 0.01;  0];
%x0_2_2 = [0.01; 0.01; 0.01; 0.01];

%a2_1 = 1.25;
%a2_2 = 1.5;

function printMatrix(Min, prec)
    fmt = ['%0.', num2str(prec), 'f  '];   % формат вывода
    for i = 1:size(Min,1)
        fprintf(fmt, Min(i,:));
        fprintf('\n');
    end
end

disp("K1:")
printMatrix(K1, 4)
disp("K2:")
printMatrix(K2, 4)

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

dt = 0.001;        % шаг интегрирования
t_values = 0:dt:14;

x0a = [0.01; 0; 0; 0];
x0dot_a = [0; 0.01; 0; 0];
x0phi = [0; 0; 0.01; 0];
x0dot_phi = [0; 0; 0; 0.01];

x0a_broken = [25; 0; 0; 0];
x0dot_a_broken = [0; 10; 0; 0];
x0phi_broken = [0; 0; 1.5; 0];
x0dot_phi_broken = [0; 0; 0; 5];

fprintf('Решение a: малые положительные значения [0.01, 0, 0, 0]\n');
[t1a, sol1a] = ode45(@(t,x) nonlinear_pendulum(t, x, K1*x, 0, M, m, l, g), t_values, x0a);

fprintf('Решение dot_a: малые положительные значения [0, 0.01, 0, 0]\n');
[t1dot_a, sol1dot_a] = ode45(@(t,x) nonlinear_pendulum(t, x, K1*x, 0, M, m, l, g), t_values, x0dot_a);

fprintf('Решение phi: малые значения [0, 0, 0.01, 0]\n');
[t1phi, sol1phi] = ode45(@(t,x) nonlinear_pendulum(t, x, K1*x, 0, M, m, l, g), t_values, x0phi);

fprintf('Решение dot_phi: очень малые значения [0, 0, 0, 0.01]\n');
[t1dot_phi, sol1dot_phi] = ode45(@(t,x) nonlinear_pendulum(t, x, K1*x, 0, M, m, l, g), t_values, x0dot_phi);

fprintf('Решение a_broken: отрицательные значения [2.25, 0, 0, 0]\n');
[t1a_broken, sol1a_broken] = ode45(@(t,x) nonlinear_pendulum(t, x, K1*x, 0, M, m, l, g), 0:dt:2, x0a_broken);

fprintf('Решение dot_a_broken: отрицательные значения [0, 2.25, 0, 0]\n');
[t1dot_a_broken, sol1dot_a_broken] = ode45(@(t,x) nonlinear_pendulum(t, x, K1*x, 0, M, m, l, g), 0:dt:2, x0dot_a_broken);

fprintf('Решение phi_broken: отрицательные значения [0, 0, 0.55, 0]\n');
[t1phi_broken, sol1phi_broken] = ode45(@(t,x) nonlinear_pendulum(t, x, K1*x, 0, M, m, l, g), 0:dt:2, x0phi_broken);

fprintf('Решение dot_phi_broken: отрицательные значения [0, 0, 0, 1.85]\n');
[t1dot_phi_broken, sol1dot_phi_broken] = ode45(@(t,x) nonlinear_pendulum(t, x, K1*x, 0, M, m, l, g), 0:dt:2, x0dot_phi_broken);




%% графики k1
figure;
plot(t1a, sol1a(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t1a, sol1a(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t1a, sol1a(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t1a, sol1a(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('Состояния при x_0 = [0.01, 0, 0, 0] и K_1');
legend('Location', 'northeast');
grid on;
saveas(gcf, 'images/fourth_part_reg_sost_min_a_k1.png');

figure;
plot(t1dot_a, sol1dot_a(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t1dot_a, sol1dot_a(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t1dot_a, sol1dot_a(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t1dot_a, sol1dot_a(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('Состояния при x_0 = [0, 0.01, 0, 0] и K_1');
legend('Location', 'northeast');
grid on;
saveas(gcf, 'images/fourth_part_reg_sost_min_dot_a_k1.png');

figure;
plot(t1phi, sol1phi(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t1phi, sol1phi(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t1phi, sol1phi(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t1phi, sol1phi(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('Состояния при x_0 = [0, 0, 0.01, 0] и K_1');
legend('Location', 'southeast');
grid on;
saveas(gcf, 'images/fourth_part_reg_sost_min_phi_k1.png');

figure;
plot(t1dot_phi, sol1dot_phi(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t1dot_phi, sol1dot_phi(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t1dot_phi, sol1dot_phi(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t1dot_phi, sol1dot_phi(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('Состояния при x_0 = [0, 0, 0, 0.01] и K_1');
legend('Location', 'northeast');
ylim([-0.025, 0.0175]);
grid on;
saveas(gcf, 'images/fourth_part_reg_sost_min_dot_phi_k1.png');

figure;
plot(t1a_broken, sol1a_broken(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t1a_broken, sol1a_broken(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t1a_broken, sol1a_broken(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t1a_broken, sol1a_broken(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('Состояния при x_0 = [25, 0, 0, 0] и K_1');
legend('Location', 'southwest');
grid on;
saveas(gcf, 'images/fourth_part_reg_sost_min_a_broken_k1.png');

figure;
plot(t1dot_a_broken, sol1dot_a_broken(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t1dot_a_broken, sol1dot_a_broken(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t1dot_a_broken, sol1dot_a_broken(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t1dot_a_broken, sol1dot_a_broken(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('Состояния при x_0 = [0, 10, 0, 0] и K_1');
legend('Location', 'southwest');
grid on;
saveas(gcf, 'images/fourth_part_reg_sost_min_dot_a_broken_k1.png');

figure;
plot(t1phi_broken, sol1phi_broken(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t1phi_broken, sol1phi_broken(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t1phi_broken, sol1phi_broken(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t1phi_broken, sol1phi_broken(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('Состояния при x_0 = [0, 0, 1.5, 0] и K_1');
legend('Location', 'southwest');
grid on;
saveas(gcf, 'images/fourth_part_reg_sost_min_phi_broken_k1.png');

figure;
plot(t1dot_phi_broken, sol1dot_phi_broken(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t1dot_phi_broken, sol1dot_phi_broken(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t1dot_phi_broken, sol1dot_phi_broken(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t1dot_phi_broken, sol1dot_phi_broken(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('Состояния при x_0 = [0, 0, 0, 5] и K_1');
legend('Location', 'northwest');
grid on;
saveas(gcf, 'images/fourth_part_reg_sost_min_dot_phi_broken_k1.png');



%% графики k2

t_values = 0:dt:6;

fprintf('Решение a: малые положительные значения [0.01, 0, 0, 0]\n');
[t2a, sol2a] = ode45(@(t,x) nonlinear_pendulum(t, x, K2*x, 0, M, m, l, g), t_values, x0a);

fprintf('Решение dot_a: малые положительные значения [0, 0.01, 0, 0]\n');
[t2dot_a, sol2dot_a] = ode45(@(t,x) nonlinear_pendulum(t, x, K2*x, 0, M, m, l, g), t_values, x0dot_a);

fprintf('Решение phi: малые значения [0, 0, 0.01, 0]\n');
[t2phi, sol2phi] = ode45(@(t,x) nonlinear_pendulum(t, x, K2*x, 0, M, m, l, g), t_values, x0phi);

fprintf('Решение dot_phi: очень малые значения [0, 0, 0, 0.01]\n');
[t2dot_phi, sol2dot_phi] = ode45(@(t,x) nonlinear_pendulum(t, x, K2*x, 0, M, m, l, g), t_values, x0dot_phi);

fprintf('Решение a_broken: отрицательные значения [2.25, 0, 0, 0]\n');
[t2a_broken, sol2a_broken] = ode45(@(t,x) nonlinear_pendulum(t, x, K2*x, 0, M, m, l, g), 0:dt:1, x0a_broken);

fprintf('Решение dot_a_broken: отрицательные значения [0, 2.25, 0, 0]\n');
[t2dot_a_broken, sol2dot_a_broken] = ode45(@(t,x) nonlinear_pendulum(t, x, K2*x, 0, M, m, l, g), 0:dt:1, x0dot_a_broken);

fprintf('Решение phi_broken: отрицательные значения [0, 0, 0.55, 0]\n');
[t2phi_broken, sol2phi_broken] = ode45(@(t,x) nonlinear_pendulum(t, x, K2*x, 0, M, m, l, g), 0:dt:1, x0phi_broken);

fprintf('Решение dot_phi_broken: отрицательные значения [0, 0, 0, 1.85]\n');
[t2dot_phi_broken, sol2dot_phi_broken] = ode45(@(t,x) nonlinear_pendulum(t, x, K2*x, 0, M, m, l, g), 0:dt:1, x0dot_phi_broken);


figure;
plot(t2a, sol2a(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t2a, sol2a(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t2a, sol2a(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t2a, sol2a(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('Состояния при x_0 = [0.01, 0, 0, 0] и K_2');
legend('Location', 'northeast');
grid on;
saveas(gcf, 'images/fourth_part_reg_sost_min_a_k2.png');

figure;
plot(t2dot_a, sol2dot_a(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t2dot_a, sol2dot_a(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t2dot_a, sol2dot_a(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t2dot_a, sol2dot_a(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('Состояния при x_0 = [0, 0.01, 0, 0] и K_2');
legend('Location', 'northeast');
grid on;
saveas(gcf, 'images/fourth_part_reg_sost_min_dot_a_k2.png');

figure;
plot(t2phi, sol2phi(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t2phi, sol2phi(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t2phi, sol2phi(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t2phi, sol2phi(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('Состояния при x_0 = [0, 0, 0.01, 0] и K_2');
legend('Location', 'southeast');
grid on;
saveas(gcf, 'images/fourth_part_reg_sost_min_phi_k2.png');

figure;
plot(t2dot_phi, sol2dot_phi(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t2dot_phi, sol2dot_phi(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t2dot_phi, sol2dot_phi(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t2dot_phi, sol2dot_phi(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('Состояния при x_0 = [0, 0, 0, 0.01] и K_2');
legend('Location', 'northeast');
ylim([-0.025, 0.0175]);
grid on;
saveas(gcf, 'images/fourth_part_reg_sost_min_dot_phi_k2.png');

figure;
plot(t2a_broken, sol2a_broken(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t2a_broken, sol2a_broken(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t2a_broken, sol2a_broken(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t2a_broken, sol2a_broken(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('Состояния при x_0 = [25, 0, 0, 0] и K_2');
legend('Location', 'southwest');
grid on;
saveas(gcf, 'images/fourth_part_reg_sost_min_a_broken_k2.png');

figure;
plot(t2dot_a_broken, sol2dot_a_broken(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t2dot_a_broken, sol2dot_a_broken(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t2dot_a_broken, sol2dot_a_broken(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t2dot_a_broken, sol2dot_a_broken(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('Состояния при x_0 = [0, 10, 0, 0] и K_2');
legend('Location', 'southwest');
grid on;
saveas(gcf, 'images/fourth_part_reg_sost_min_dot_a_broken_k2.png');

figure;
plot(t2phi_broken, sol2phi_broken(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t2phi_broken, sol2phi_broken(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t2phi_broken, sol2phi_broken(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t2phi_broken, sol2phi_broken(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('Состояния при x_0 = [0, 0, 1.5, 0] и K_2');
legend('Location', 'southwest');
grid on;
saveas(gcf, 'images/fourth_part_reg_sost_min_phi_broken_k2.png');

figure;
plot(t2dot_phi_broken, sol2dot_phi_broken(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t2dot_phi_broken, sol2dot_phi_broken(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t2dot_phi_broken, sol2dot_phi_broken(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t2dot_phi_broken, sol2dot_phi_broken(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('Состояния при x_0 = [0, 0, 0, 5] и K_2');
legend('Location', 'northwest');
grid on;
saveas(gcf, 'images/fourth_part_reg_sost_min_dot_phi_broken_k2.png');