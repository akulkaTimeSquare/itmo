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

a1 = 0.5;
a2 = 1.5;
a3 = 3;
a4 = 6;

cvx_begin sdp
variable P1(4, 4) symmetric
variable Y1(1, 4)
P1 > 0.0001*eye(4);
P1*A' + A*P1 + 2*a1*P1 + Y1'*B' + B*Y1 <= 0;

variable P2(4, 4) symmetric
variable Y2(1, 4)
P2 > 0.0001*eye(4);
P2*A' + A*P2 + 2*a2*P2 + Y2'*B' + B*Y2 <= 0;

variable P3(4, 4) symmetric
variable Y3(1, 4)
P3 > 0.0001*eye(4);
P3*A' + A*P3 + 2*a3*P3 + Y3'*B' + B*Y3 <= 0;

variable P4(4, 4) symmetric
variable Y4(1, 4)
P4 > 0.0001*eye(4);
P4*A' + A*P4 + 2*a4*P4 + Y4'*B' + B*Y4 <= 0;

cvx_end

K1 = Y1/P1;
K2 = Y2/P2;
K3 = Y3/P3;
K4 = Y4/P4;

e1 = eig(A + B*K1)
e2 = eig(A + B*K2)
e3 = eig(A + B*K3)
e4 = eig(A + B*K4)

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
disp("K3:")
printMatrix(K3, 4)
disp("K4:")
printMatrix(K4, 4)

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

t_span = [0 5];  % временной интервал
dt = 0.0001;        % шаг интегрирования
t_values = 0:dt:35;

% Начальные условия из CS_project (все состояния ненулевые)
x0 = [-0.1; 0.05; -0.15; 0.075];

fprintf('Решение 1');
[t1, sol1] = ode45(@(t,x) nonlinear_pendulum(t, x, K1*x, 0, M, m, l, g), t_values, x0);
u1 = K1*sol1';
phi1max = max(abs(sol1(:,3)))
a1 = max(abs(sol1(:,1)))
u1max = max(abs(u1))

fprintf('Решение 2');
[t2, sol2] = ode45(@(t,x) nonlinear_pendulum(t, x, K2*x, 0, M, m, l, g), t_values, x0);
u2 = K2*sol2';
phi2max = max(abs(sol2(:,3)))
a2 = max(abs(sol2(:,1)))
u2max = max(abs(u2))

fprintf('Решение 3');
[t3, sol3] = ode45(@(t,x) nonlinear_pendulum(t, x, K3*x, 0, M, m, l, g), t_values, x0);
u3 = K3*sol3';
phi3max = max(abs(sol3(:,3)))
a3 = max(abs(sol3(:,1)))
u3max = max(abs(u3))

fprintf('Решение 4');
[t4, sol4] = ode45(@(t,x) nonlinear_pendulum(t, x, K4*x, 0, M, m, l, g), t_values, x0);
u4 = K4*sol4';
phi4max = max(abs(sol4(:,3)))
a4 = max(abs(sol4(:,1)))
u4max = max(abs(u4))

figure;
plot(t1, sol1(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t1, sol1(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t1, sol1(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t1, sol1(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('x(t) нелинейной модели при \alpha_1 = 0.5');
legend('Location', 'northeast');
grid on;
xlim([-5e-3, 10])
saveas(gcf, 'images/fourth_part_reg_sost_an_1.png');

figure;
plot(t1, u1, 'LineWidth', 2, "Clipping", "on");
xlabel('t');
ylabel('u');
title('Управляющее воздействие при \alpha_1 = 0.5');
grid on;
xlim([-5e-3, 10])
ylim([-1200, 6000])
saveas(gcf, 'images/fourth_part_reg_sost_an_1_u.png');

figure;
plot(t2, sol2(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t2, sol2(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t2, sol2(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t2, sol2(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('x(t) нелинейной модели при \alpha_2 = 1.5');
legend('Location', 'northeast');
grid on;
xlim([-5e-3, 3])
saveas(gcf, 'images/fourth_part_reg_sost_an_2.png');

figure;
plot(t2, u2, 'LineWidth', 2, "Clipping", "on");
xlabel('t');
ylabel('u');
title('Управляющее воздействие при \alpha_2 = 1.5');
grid on;
xlim([-5e-3, 3])
saveas(gcf, 'images/fourth_part_reg_sost_an_2_u.png');

figure;
plot(t3, sol3(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t3, sol3(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t3, sol3(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t3, sol3(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('x(t) нелинейной модели при \alpha_3 = 3');
legend('Location', 'northeast');
grid on;
xlim([-5e-3, 3])
ylim([-0.6, 2.25])
saveas(gcf, 'images/fourth_part_reg_sost_an_3.png');

figure;
plot(t3, u3, 'LineWidth', 2, "Clipping", "on");
xlabel('t');
ylabel('u');
title('Управляющее воздействие при \alpha_3 = 3');
grid on;
xlim([-5e-3, 3])
ylim([-1.25*10^4, 5*10^4])
saveas(gcf, 'images/fourth_part_reg_sost_an_3_u.png');

figure;
plot(t4, sol4(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t4, sol4(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t4, sol4(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t4, sol4(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('x(t) нелинейной модели при \alpha_4 = 6');
legend('Location', 'northeast');
grid on;
xlim([-5e-3, 1])
saveas(gcf, 'images/fourth_part_reg_sost_an_4.png');

figure;
plot(t4, u4, 'LineWidth', 2, "Clipping", "on");
xlabel('t');
ylabel('u');
title('Управляющее воздействие при \alpha_4 = 6');
grid on;
xlim([-5e-3, 1])
saveas(gcf, 'images/fourth_part_reg_sost_an_4_u.png');