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

G1 = [-0.25  0 0 0; 0 -0.5 0 0; 0 0 -0.75 0; 0 0 0 -1];
G2 = [-2 1 0 0; 0 -2 0 0; 0 0 -3 1; 0 0  0 -3];
G3 = [-2 4 0 0; -4 -2 0 0; 0 0 -3 4; 0 0 -4 -3];
G4 = [-4  0 0 0; 0 -6 0 0; 0 0 -8 0; 0 0 0 -10];
Y = [1 1 1 1];

rank([Y; Y*G1; Y*G1^2; Y*G1^3])
rank([Y; Y*G2; Y*G2^2; Y*G2^3])
rank([Y; Y*G3; Y*G3^2; Y*G3^3])
rank([Y; Y*G4; Y*G4^2; Y*G4^3])

P1 = sylvester(A, -G1, B*Y);
P2 = sylvester(A, -G2, B*Y);
P3 = sylvester(A, -G3, B*Y);
P4 = sylvester(A, -G4, B*Y);

K1 = -Y*P1^(-1);
K2 = -Y*P2^(-1);
K3 = -Y*P3^(-1);
K4 = -Y*P4^(-1);

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

disp("K1:");
printMatrix(K1, 4);
disp("K2:");
printMatrix(K2, 4);
disp("K3:");
printMatrix(K3, 4);
disp("K4:");
printMatrix(K4, 4);

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
x0 = [-0.01; -0.2; 0.3; -0.15];

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

%% графики'
figure;
plot(t1, sol1(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t1, sol1(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t1, sol1(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t1, sol1(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('x(t) нелинейной модели при K_1');
legend('Location', 'southeast');
grid on;
xlim([0, 35])
saveas(gcf, 'images/third_part_reg_sost_an_1.png');

figure;
plot(t1, u1, 'LineWidth', 2, "Clipping", "on");
xlabel('t');
ylabel('u');
title('Управляющее воздействие при K_1');
grid on;
xlim([0, 35])
saveas(gcf, 'images/third_part_reg_sost_an_1_u.png');

figure;
plot(t2, sol2(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t2, sol2(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t2, sol2(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t2, sol2(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('x(t) нелинейной модели при K_2');
legend('Location', 'southeast');
grid on;
xlim([0, 5])
saveas(gcf, 'images/third_part_reg_sost_an_2.png');

figure;
plot(t2, u2, 'LineWidth', 2, "Clipping", "on");
xlabel('t');
ylabel('u');
title('Управляющее воздействие при K_2');
grid on;
xlim([0, 5])
saveas(gcf, 'images/third_part_reg_sost_an_2_u.png');

figure;
plot(t3, sol3(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t3, sol3(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t3, sol3(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t3, sol3(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('x(t) нелинейной модели при K_3');
legend('Location', 'northeast');
grid on;
xlim([0, 5])
saveas(gcf, 'images/third_part_reg_sost_an_3.png');

figure;
plot(t3, u3, 'LineWidth', 2, "Clipping", "on");
xlabel('t');
ylabel('u');
title('Управляющее воздействие при K_3');
grid on;
xlim([0, 5])
saveas(gcf, 'images/third_part_reg_sost_an_3_u.png');

figure;
plot(t4, sol4(:,1), 'LineWidth', 2, 'DisplayName', 'x_1', "Clipping", "on"); hold on;
plot(t4, sol4(:,2), 'LineWidth', 2, 'DisplayName', 'x_2', "Clipping", "on"); 
plot(t4, sol4(:,3), 'LineWidth', 2, 'DisplayName', 'x_3', "Clipping", "on"); 
plot(t4, sol4(:,4), 'LineWidth', 2, 'DisplayName', 'x_4', "Clipping", "on");
xlabel('t');
ylabel('x');
title('x(t) нелинейной модели при K_4');
legend('Location', 'northeast');
ylim([-6, 3.25])
xlim([0, 3])
grid on;
saveas(gcf, 'images/third_part_reg_sost_an_4.png');

figure;
plot(t4, u4, 'LineWidth', 2, "Clipping", "on");
xlabel('t');
ylabel('u');
title('Управляющее воздействие при K_4');
grid on;
xlim([0, 3])
saveas(gcf, 'images/third_part_reg_sost_an_4_u.png');