%%
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

G = [-0.25 0 0 0; 0 -0.5 0 0; 0 0 -0.75 0; 0 0  0 -1];
Yr = [1 1 1 1];

rank([Yr; Yr*G; Yr*G^2; Yr*G^3])

P = sylvester(A, -G, B*Yr);
K = -Yr*P^(-1);
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

Gl1 = [-0.75  0 0 0; 0 -1.25 0 0; 0 0 -1.5 0; 0 0 0 -2];
Gl2 = [-2 1 0 0; 0 -2 0 0; 0 0 -4 1; 0 0  0 -4];
Gl3 = [-2 4 0 0; -4 -2 0 0; 0 0 -4 4; 0 0 -4 -4];
Gl4 = [-7  0 0 0; 0 -7.5 0 0; 0 0 -8 0; 0 0 0 -8.5];

Yl = [1 1; 1 1; 1 1; 1 1];
V1 = [Yl Gl1*Yl Gl1^2*Yl Gl1^3*Yl];
V2 = [Yl Gl2*Yl Gl2^2*Yl Gl2^3*Yl];
V3 = [Yl Gl3*Yl Gl3^2*Yl Gl3^3*Yl];
V4 = [Yl Gl4*Yl Gl4^2*Yl Gl4^3*Yl];
display(rank(V1))
display(rank(V2))
display(rank(V3))
display(rank(V4))

Q1 = sylvester(Gl1, -A, Yl*C);
L1 = inv(Q1)*Yl
Q2 = sylvester(Gl2, -A, Yl*C);
L2 = inv(Q2)*Yl
Q3 = sylvester(Gl3, -A, Yl*C);
L3 = inv(Q3)*Yl
Q4 = sylvester(Gl4, -A, Yl*C);
L4 = inv(Q4)*Yl

e1 = eig(A + L1*C)
e2 = eig(A + L2*C)
e3 = eig(A + L3*C)
e4 = eig(A + L4*C)

disp("L1:");
printMatrix(L1, 4);
disp("L2:");
printMatrix(L2, 4);
disp("L3:");
printMatrix(L3, 4);
disp("L4:");
printMatrix(L4, 4);

% Nonlinear system simulation with observer
dt = 0.0001;
t_values = 0:dt:25;

% Define the nonlinear system function
function dzdt = nonlinear_system_observer(t, z, L, A, B, C, K, M, m, l, g)
    x = z(1:4);
    x_hat = z(5:8);
    
    y_hat = C * x_hat;
    u = (K * x);
    y = C * x;
    f_val = 0;
    
    % Nonlinear dynamics
    denominator1 = (1/3)*(M + m)*l - (1/4)*m*l*cos(x(3))^2;
    numerator1 = (1/3)*l*(u - 0.5*m*l*x(4)^2*sin(x(3))) + 0.5*cos(x(3))*(f_val + 0.5*m*g*l*sin(x(3)));
    dx2 = numerator1 / denominator1;
    
    denominator2 = (1/3)*(M + m)*m*l^2 - (1/4)*m^2*l^2*cos(x(3))^2;
    numerator2 = 0.5*m*l*cos(x(3))*(u - 0.5*m*l*x(4)^2*sin(x(3))) + (M + m)*(f_val + 0.5*m*g*l*sin(x(3)));
    dx4 = numerator2 / denominator2;
    
    % Observer dynamics
    dx_hat = A * x_hat + B * u + L * (y_hat - y);
    
    dzdt = [x(2); dx2; x(4); dx4; dx_hat];
end

% Initial conditions
x0 = [0.025; 0.03; 0.015; -0.01];
x_hat0 = [0; 0; 0; 0];
z0 = [x0; x_hat0];

% Solve for different observers
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

[t1, sol_1] = ode45(@(t,z) nonlinear_system_observer(t, z, L1, A, B, C, K, M, m, l, g), t_values, z0, options);
[t2, sol_2] = ode45(@(t,z) nonlinear_system_observer(t, z, L2, A, B, C, K, M, m, l, g), t_values, z0, options);
[t3, sol_3] = ode45(@(t,z) nonlinear_system_observer(t, z, L3, A, B, C, K, M, m, l, g), t_values, z0, options);
[t4, sol_4] = ode45(@(t,z) nonlinear_system_observer(t, z, L4, A, B, C, K, M, m, l, g), t_values, z0, options);

% Extract real and estimated states
x_real_1 = sol_1(:, 1:4);
x_real_2 = sol_2(:, 1:4);
x_real_3 = sol_3(:, 1:4);
x_real_4 = sol_4(:, 1:4);

x_est_1 = sol_1(:, 5:8);
x_est_2 = sol_2(:, 5:8);
x_est_3 = sol_3(:, 5:8);
x_est_4 = sol_4(:, 5:8);

% Calculate errors
e_1 = x_real_1 - x_est_1;
e_2 = x_real_2 - x_est_2;
e_3 = x_real_3 - x_est_3;
e_4 = x_real_4 - x_est_4;

figure;
plot(t1, e_1(:,1), "LineWidth", 2);
hold on;
plot(t1, e_1(:,2), "LineWidth", 2);
plot(t1, e_1(:,3), "LineWidth", 2);
plot(t1, e_1(:,4), "LineWidth", 2);
title("Ошибка оценки e_1 при L_1");
grid on;
xlabel("t");
ylabel("e_1");
xlim([0, 15]);
legend("e_1(1)", "e_1(2)", "e_1(3)", "e_1(4)", "Location", "southeast");
saveas(gcf, "images/third_part_observer_e_1.png");

figure;
plot(t1, x_real_1(:,1), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t1, x_est_1(:,1), "LineWidth", 3, LineStyle="-.");
title("x_1(t) при L_1", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_1", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "southeast");
set(gca, "FontSize", 14);
saveas(gcf, "images/third_part_observer1_x_1.png");

figure;
plot(t1, x_real_1(:,2), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t1, x_est_1(:,2), "LineWidth", 3, LineStyle="-.");
title("x_2(t) при L_1", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_2", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13);
set(gca, "FontSize", 14);
ylim([-0.105, 0.04]);
saveas(gcf, "images/third_part_observer1_x_2.png");

figure;
plot(t1, x_real_1(:,3), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t1, x_est_1(:,3), "LineWidth", 3, LineStyle="-.");
title("x_3(t) при L_1", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_3", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13);
set(gca, "FontSize", 14);
xlim([0, 12]);
saveas(gcf, "images/third_part_observer1_x_3.png");

figure;
plot(t1, x_real_1(:,4), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t1, x_est_1(:,4), "LineWidth", 3, LineStyle="-.");
title("x_4(t) при L_1", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_4", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13);
set(gca, "FontSize", 14);
xlim([0, 12]);
saveas(gcf, "images/third_part_observer1_x_4.png");

figure;
plot(t2, e_2(:,1), "LineWidth", 2);
hold on;
plot(t2, e_2(:,2), "LineWidth", 2);
plot(t2, e_2(:,3), "LineWidth", 2);
plot(t2, e_2(:,4), "LineWidth", 2);
title("Ошибка оценки e_2 при L_2");
grid on;
xlabel("t");
ylabel("e_2");
legend("e_2(1)", "e_2(2)", "e_2(3)", "e_2(4)", "Location", "southeast");
xlim([0, 5]);
saveas(gcf, "images/third_part_observer_e_2.png");

figure;
plot(t2, x_real_2(:,1), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t2, x_est_2(:,1), "LineWidth", 3, LineStyle="-.");
title("x_1(t) при L_2", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_1", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "northeast");
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/third_part_observer2_x_1.png");

figure;
plot(t2, x_real_2(:,2), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t2, x_est_2(:,2), "LineWidth", 3, LineStyle="-.");
title("x_2(t) при L_2", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_2", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "southeast");
set(gca, "FontSize", 14);
ylim([-0.082, 0.04]);
xlim([0, 5]);
saveas(gcf, "images/third_part_observer2_x_2.png");

figure;
plot(t2, x_real_2(:,3), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t2, x_est_2(:,3), "LineWidth", 3, LineStyle="-.");
title("x_3(t) при L_2", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_3", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13);
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/third_part_observer2_x_3.png");

figure;
plot(t2, x_real_2(:,4), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t2, x_est_2(:,4), "LineWidth", 3, LineStyle="-.");
title("x_4(t) при L_2", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_4", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13);
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/third_part_observer2_x_4.png");

figure;
plot(t3, e_3(:,1), "LineWidth", 2);
hold on;
plot(t3, e_3(:,2), "LineWidth", 2);
plot(t3, e_3(:,3), "LineWidth", 2);
plot(t3, e_3(:,4), "LineWidth", 2);
title("Ошибка оценки e_3 при L_3");
grid on;
xlabel("t");
ylabel("e_3");
xlim([0, 5]);
legend("e_3(1)", "e_3(2)", "e_3(3)", "e_3(4)", "Location", "southeast");
saveas(gcf, "images/third_part_observer_e_3.png");

figure;
plot(t3, x_real_3(:,1), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t3, x_est_3(:,1), "LineWidth", 3, LineStyle="-.");
title("x_1(t) при L_3", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_1", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13);
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/third_part_observer3_x_1.png");

figure;
plot(t3, x_real_3(:,2), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t3, x_est_3(:,2), "LineWidth", 3, LineStyle="-.");
title("x_2(t) при L_3", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_2", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13);
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/third_part_observer3_x_2.png");

figure;
plot(t3, x_real_3(:,3), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t3, x_est_3(:,3), "LineWidth", 3, LineStyle="-.");
title("x_3(t) при L_3", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_3", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13);
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/third_part_observer3_x_3.png");

figure;
plot(t3, x_real_3(:,4), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t3, x_est_3(:,4), "LineWidth", 3, LineStyle="-.");
title("x_4(t) при L_3", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_4", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13);
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/third_part_observer3_x_4.png");

figure;
plot(t4, e_4(:,1), "LineWidth", 2);
hold on;
plot(t4, e_4(:,2), "LineWidth", 2);
plot(t4, e_4(:,3), "LineWidth", 2);
plot(t4, e_4(:,4), "LineWidth", 2);
title("Ошибка оценки e_4 при L_4");
grid on;
xlabel("t");
ylabel("e_4");
xlim([0, 3]);
legend("e_4(1)", "e_4(2)", "e_4(3)", "e_4(4)", "Location", "southeast");
saveas(gcf, "images/third_part_observer_e_4.png");

%%
figure;
plot(t4, x_real_4(:,1), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t4, x_est_4(:,1), "LineWidth", 3, LineStyle="-.");
title("x_1(t) при L_4", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_1", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13);
set(gca, "FontSize", 14);
xlim([0, 3]);
saveas(gcf, "images/third_part_observer4_x_1.png");

%%
figure;
plot(t4, x_real_4(:,2), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t4, x_est_4(:,2), "LineWidth", 3, LineStyle="-.");
title("x_2(t) при L_4", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_2", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13);
set(gca, "FontSize", 14);
xlim([0, 3]);
saveas(gcf, "images/third_part_observer4_x_2.png");

figure;
plot(t4, x_real_4(:,3), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t4, x_est_4(:,3), "LineWidth", 3, LineStyle="-.");
title("x_3(t) при L_4", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_3", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13);
set(gca, "FontSize", 14);
xlim([0, 3]);
saveas(gcf, "images/third_part_observer4_x_3.png");

figure;
plot(t4, x_real_4(:,4), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t4, x_est_4(:,4), "LineWidth", 3, LineStyle="-.");
title("x_4(t) при L_4", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_4", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13);
set(gca, "FontSize", 14);
xlim([0, 3]);
saveas(gcf, "images/third_part_observer4_x_4.png");

%%
% Reduced-order observer simulation
dt = 0.0001;
x0 = [0.025; 0.03; 0.015; -0.01];
t_eval = 0:dt:5;

% Define the reduced-order nonlinear system function
function dzdt = nonlinear_system_reduced(t, z, Q, G, Y, A, B, C, K, M, m, l, g)
    x = z(1:4);
    z_hat = z(5:6);  % Reduced order observer state (2D)

    y = C * x;
    u = K * x;

    f_val = 0;
    
    % Nonlinear dynamics
    denominator1 = (1/3)*(M + m)*l - (1/4)*m*l*cos(x(3))^2;
    numerator1 = (1/3)*l*(u - 0.5*m*l*x(4)^2*sin(x(3))) + 0.5*cos(x(3))*(f_val + 0.5*m*g*l*sin(x(3)));
    dx2 = numerator1 / denominator1;
    
    denominator2 = (1/3)*(M + m)*m*l^2 - (1/4)*m^2*l^2*cos(x(3))^2;
    numerator2 = 0.5*m*l*cos(x(3))*(u - 0.5*m*l*x(4)^2*sin(x(3))) + (M + m)*(f_val + 0.5*m*g*l*sin(x(3)));
    dx4 = numerator2 / denominator2;
    
    % Reduced observer dynamics
    dzhat_dt = G * z_hat - Y * y + (Q * (B * u));
    
    dzdt = [x(2); dx2; x(4); dx4; dzhat_dt];
end

% Initial conditions for reduced observer
x0_red = x0;
x_hat0_red = [0; 0];  % 2D observer state
z0_red = [x0_red; x_hat0_red];

% Solve for reduced observers (using the 2D observers from earlier)
options_red = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% We need to define Q1, Q2, G1, G2, Y1, Y2 for the reduced observers
Gd1 = [-0.5 0; 0 -0.75];
Gd2 = [-3 0; 0 -4];

Yd = [-1 1; 1 0];
Ud1 = [Yd, Gd1*Yd];
Ud2 = [Yd, Gd2*Yd];
rankUd1 = rank(Ud1)
rankUd2 = rank(Ud2)

Qd1 = sylvester(Gd1, -A, Yd*C);
Qd2 = sylvester(Gd2, -A, Yd*C);

disp("Qd1:");
disp(Qd1);
disp("Qd2:");
disp(Qd2);

CQd1 = [C; Qd1];
CQd2 = [C; Qd2];
invCQd1 = inv(CQd1);
invCQd2 = inv(CQd2);

[t_red1, sol_red1] = ode45(@(t,z) nonlinear_system_reduced(t, z, Qd1, Gd1, Yd, A, B, C, K, M, m, l, g), 0:dt:30, z0_red, options_red);
[t_red2, sol_red2] = ode45(@(t,z) nonlinear_system_reduced(t, z, Qd2, Gd2, Yd, A, B, C, K, M, m, l, g), 0:dt:30, z0_red, options_red);

% Extract solutions
x_sol_1_red = sol_red1(:, 1:4)';  % Transpose to match Python indexing
z_sol_1_red = sol_red1(:, 5:6)';
x_sol_2_red = sol_red2(:, 1:4)';
z_sol_2_red = sol_red2(:, 5:6)';

y_sol_1_red = C * x_sol_1_red;
y_sol_2_red = C * x_sol_2_red;

% Reconstruct full state estimates
x_hat_sol_1_red = invCQd1 * [y_sol_1_red; z_sol_1_red];
x_hat_sol_2_red = invCQd2 * [y_sol_2_red; z_sol_2_red];

% Calculate errors for reduced observers
e_sol_1_red = x_sol_1_red - x_hat_sol_1_red;
e_sol_2_red = x_sol_2_red - x_hat_sol_2_red;

figure;
plot(t_red1, e_sol_1_red(1,:), "LineWidth", 2);
hold on;
plot(t_red1, e_sol_1_red(2,:), "LineWidth", 2);
plot(t_red1, e_sol_1_red(3,:), "LineWidth", 2);
plot(t_red1, e_sol_1_red(4,:), "LineWidth", 2);
title("Ошибка оценки e_1 при Q_1");
grid on;
xlabel("t");
ylabel("e_1");
xlim([0, 15]);
ylim([-0.025, 0.2])
legend("e_1(1)", "e_1(2)", "e_1(3)", "e_1(4)", "Location", "northeast");
saveas(gcf, "images/third_part_observer_red_e_1.png");

figure;
plot(t_red1, x_sol_1_red(1,:), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t_red1, x_hat_sol_1_red(1,:), "LineWidth", 3, LineStyle="-.");
title("x_1(t) при Q_1", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_1", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "southeast");
set(gca, "FontSize", 14);
saveas(gcf, "images/third_part_observer1_red_x_1.png");

figure;
plot(t_red1, x_sol_1_red(2,:), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t_red1, x_hat_sol_1_red(2,:), "LineWidth", 3, LineStyle="-.");
title("x_2(t) при Q_1", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_2", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13);
set(gca, "FontSize", 14);
saveas(gcf, "images/third_part_observer1_red_x_2.png");

figure;
plot(t_red1, x_sol_1_red(3,:), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t_red1, x_hat_sol_1_red(3,:), "LineWidth", 3, LineStyle="-.");
title("x_3(t) при Q_1", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_3", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13);
set(gca, "FontSize", 14);
saveas(gcf, "images/third_part_observer1_red_x_3.png");

figure;
plot(t_red1, x_sol_1_red(4,:), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t_red1, x_hat_sol_1_red(4,:), "LineWidth", 3, LineStyle="-.");
title("x_4(t) при Q_1", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_4", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "southeast");
set(gca, "FontSize", 14);
saveas(gcf, "images/third_part_observer1_red_x_4.png");


%%
figure;
plot(t_red2, e_sol_2_red(1,:), "LineWidth", 2);
hold on;
plot(t_red2, e_sol_2_red(2,:), "LineWidth", 2);
plot(t_red2, e_sol_2_red(3,:), "LineWidth", 2);
plot(t_red2, e_sol_2_red(4,:), "LineWidth", 2);
title("Ошибка оценки e_2 при Q_2");
grid on;
xlabel("t");
ylabel("e_1");
xlim([0, 5]);
ylim([-0.08, 0.01]);
legend("e_2(1)", "e_2(2)", "e_2(3)", "e_2(4)", "Location", "southeast");
saveas(gcf, "images/third_part_observer_red_e_2.png");

figure;
plot(t_red2, x_sol_2_red(1,:), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t_red2, x_hat_sol_2_red(1,:), "LineWidth", 3, LineStyle="-.");
title("x_1(t) при Q_2", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_1", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "northeast");
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/third_part_observer2_red_x_1.png");

figure;
plot(t_red2, x_sol_2_red(2,:), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t_red2, x_hat_sol_2_red(2,:), "LineWidth", 3, LineStyle="-.");
title("x_2(t) при Q_2", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_2", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13);
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/third_part_observer2_red_x_2.png");

figure;
plot(t_red2, x_sol_2_red(3,:), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t_red2, x_hat_sol_2_red(3,:), "LineWidth", 3, LineStyle="-.");
title("x_3(t) при Q_2", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_3", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13);
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/third_part_observer2_red_x_3.png");

figure;
plot(t_red2, x_sol_2_red(4,:), "LineWidth", 3, "Clipping", "on"); hold on;
plot(t_red2, x_hat_sol_2_red(4,:), "LineWidth", 3, LineStyle="-.");
title("x_4(t) при Q_2", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_4", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13);
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/third_part_observer2_red_x_4.png");