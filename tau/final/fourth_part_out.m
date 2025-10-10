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


ar1 = 1;
ar2 = 1;
ar3 = 4;
ar4 = 4;

al1 = 1;
al2 = 4;
al3 = 1;
al4 = 4;

cvx_begin sdp quiet
    variable Pr1(4, 4) symmetric
    variable Yr1(1, 4)
    Pr1 > 0.0001*eye(4);
    Pr1*A' + A*Pr1 + 2*ar1*Pr1 + Yr1'*B' + B*Yr1 <= 0;

    variable Ql1(4, 4) symmetric
    variable Yl1(4, 2)
    Ql1 > 0.0001*eye(4);
    A'*Ql1 + Ql1*A + 2*al1*Ql1 + C'*Yl1' + Yl1*C <= 0;

    variable Pr2(4, 4) symmetric
    variable Yr2(1, 4)
    Pr2 > 0.0001*eye(4);
    Pr2*A' + A*Pr2 + 2*ar2*Pr2 + Yr2'*B' + B*Yr2 <= 0;

    variable Ql2(4, 4) symmetric
    variable Yl2(4, 2)
    Ql2 > 0.0001*eye(4);
    A'*Ql2 + Ql2*A + 2*al2*Ql2 + C'*Yl2' + Yl2*C <= 0;

    variable Pr3(4, 4) symmetric
    variable Yr3(1, 4)
    Pr3 > 0.0001*eye(4);
    Pr3*A' + A*Pr3 + 2*ar3*Pr3 + Yr3'*B' + B*Yr3 <= 0;

    variable Ql3(4, 4) symmetric
    variable Yl3(4, 2)
    Ql3 > 0.0001*eye(4);
    A'*Ql3 + Ql3*A + 2*al3*Ql3 + C'*Yl3' + Yl3*C <= 0;

    variable Pr4(4, 4) symmetric
    variable Yr4(1, 4)
    Pr4 > 0.0001*eye(4);
    Pr4*A' + A*Pr4 + 2*ar4*Pr4 + Yr4'*B' + B*Yr4 <= 0;

    variable Ql4(4, 4) symmetric
    variable Yl4(4, 2)
    Ql4 > 0.0001*eye(4);
    A'*Ql4 + Ql4*A + 2*al4*Ql4 + C'*Yl4' + Yl4*C <= 0;
cvx_end

K1 = Yr1/Pr1;
K2 = Yr2/Pr2;
K3 = Yr3/Pr3;
K4 = Yr4/Pr4;

L1 = Ql1\Yl1;
L2 = Ql2\Yl2;
L3 = Ql3\Yl3;
L4 = Ql4\Yl4;

er1 = eig(A + B*K1)
er2 = eig(A + B*K2)
er3 = eig(A + B*K3)
er4 = eig(A + B*K4)

el1 = eig(A + L1*C)
el2 = eig(A + L2*C)
el3 = eig(A + L3*C)
el4 = eig(A + L4*C)

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

disp("L1:");
printMatrix(L1', 4);
disp("L2:");
printMatrix(L2', 4);
disp("L3:");
printMatrix(L3', 4);
disp("L4:");
printMatrix(L4', 4);

function dzdt = nonlinear_system_observer(t, z, L, A, B, C, K, M, m, l, g)
    x = z(1:4);
    x_hat = z(5:8);
    
    y_hat = C * x_hat;
    u = K * x_hat;
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
x0 = [0.05; 0.06; 0.07; 0.02];
x_hat0 = [0; 0; 0; 0];
z0 = [x0; x_hat0];

% Solve for different observers
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
dt = 0.0001;
t_values = 0:dt:25;

[t1, sol_1] = ode45(@(t,z) nonlinear_system_observer(t, z, L1, A, B, C, K1, M, m, l, g), t_values, z0, options);
[t2, sol_2] = ode45(@(t,z) nonlinear_system_observer(t, z, L2, A, B, C, K2, M, m, l, g), t_values, z0, options);
[t3, sol_3] = ode45(@(t,z) nonlinear_system_observer(t, z, L3, A, B, C, K3, M, m, l, g), t_values, z0, options);
[t4, sol_4] = ode45(@(t,z) nonlinear_system_observer(t, z, L4, A, B, C, K4, M, m, l, g), t_values, z0, options);

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
e_1 = (x_real_1 - x_est_1)';
e_2 = (x_real_2 - x_est_2)';
e_3 = (x_real_3 - x_est_3)';
e_4 = (x_real_4 - x_est_4)';

u1 = K1*x_est_1';
u2 = K2*x_est_2';
u3 = K3*x_est_3';
u4 = K4*x_est_4';

%%
figure;
plot(t1, u1, "LineWidth", 3);
title("Управление u_1(t) при K_1 и L_1", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("u_1", "FontSize", 15);
set(gca, "FontSize", 14);
xlim([0, 6]);
saveas(gcf, "images/fourth_part_out_u1.png");

figure;
plot(t1, e_1(1,:), "LineWidth", 3);
hold on;
plot(t1, e_1(2, :), "LineWidth", 3);
plot(t1, e_1(3, :), "LineWidth", 3);
plot(t1, e_1(4, :), "LineWidth", 3);
title("Ошибка оценки e_1(t) при K_1 и L_1", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("e_1", "FontSize", 15);
legend("e_1(1)", "e_1(2)", "e_1(3)", "e_1(4)", "FontSize", 13, "Location", "northeast");
set(gca, "FontSize", 14);
xlim([0, 6]);
saveas(gcf, "images/fourth_part_out_e_1.png");

figure;
plot(t1, x_real_1(:, 1), "LineWidth", 3);
hold on;
plot(t1, x_est_1(:, 1), "LineWidth", 3, "LineStyle", "-.");
title("x_1(t) при K_1 и L_1", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_1", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "southeast");
set(gca, "FontSize", 14);
xlim([0, 6]);
saveas(gcf, "images/fourth_part_out1_x_1.png");

figure;
plot(t1, x_real_1(:, 2), "LineWidth", 3);
hold on;
plot(t1, x_est_1(:, 2), "LineWidth", 3, "LineStyle", "-.");
title("x_2(t) при K_1 и L_1", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_2", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "northeast");
set(gca, "FontSize", 14);
xlim([0, 6]);
saveas(gcf, "images/fourth_part_out1_x_2.png");

figure;
plot(t1, x_real_1(:, 3), "LineWidth", 3);
hold on;
plot(t1, x_est_1(:, 3), "LineWidth", 3, "LineStyle", "-.");
title("x_3(t) при K_1 и L_1", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_3", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "northeast");
set(gca, "FontSize", 14);
xlim([0, 6]);
saveas(gcf, "images/fourth_part_out1_x_3.png");

figure;
plot(t1, x_real_1(:, 4), "LineWidth", 3);
hold on;
plot(t1, x_est_1(:, 4), "LineWidth", 3, "LineStyle", "-.");
title("x_4(t) при K_1 и L_1", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_4", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "northeast");
set(gca, "FontSize", 14);
xlim([0, 6]);
saveas(gcf, "images/fourth_part_out1_x_4.png");



figure;
plot(t2, u2, "LineWidth", 3);
title("Управление u_2(t) при K_2 и L_2", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("u_2", "FontSize", 15);
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/fourth_part_out2_u2.png");

figure;
plot(t2, e_2(1,:), "LineWidth", 3);
hold on;
plot(t2, e_2(2, :), "LineWidth", 3);
plot(t2, e_2(3, :), "LineWidth", 3);
plot(t2, e_2(4, :), "LineWidth", 3);
title("Ошибка оценки e_2(t) при K_2 и L_2", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("e_2", "FontSize", 15);
legend("e_2(1)", "e_2(2)", "e_2(3)", "e_2(4)", "FontSize", 13, "Location", "southeast");
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/fourth_part_out2_e_2.png");

figure;
plot(t2, x_real_2(:, 1), "LineWidth", 3);
hold on;
plot(t2, x_est_2(:, 1), "LineWidth", 3, "LineStyle", "-.");
title("x_1(t) при K_2 и L_2", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_1", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "southeast");
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/fourth_part_out2_x_1.png");

figure;
plot(t2, x_real_2(:, 2), "LineWidth", 3);
hold on;
plot(t2, x_est_2(:, 2), "LineWidth", 3, "LineStyle", "-.");
title("x_2(t) при K_2 и L_2", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_2", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "northeast");
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/fourth_part_out2_x_2.png");

figure;
plot(t2, x_real_2(:, 3), "LineWidth", 3);
hold on;
plot(t2, x_est_2(:, 3), "LineWidth", 3, "LineStyle", "-.");
title("x_3(t) при K_2 и L_2", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_3", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "northeast");
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/fourth_part_out2_x_3.png");

figure;
plot(t2, x_real_2(:, 4), "LineWidth", 3);
hold on;
plot(t2, x_est_2(:, 4), "LineWidth", 3, "LineStyle", "-.");
title("x_4(t) при K_2 и L_2", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_4", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "northeast");
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/fourth_part_out2_x_4.png");

%%
figure;
plot(t3, u3, "LineWidth", 3);
title("Управление u_3(t) при K_3 и L_3", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("u_3", "FontSize", 15);
set(gca, "FontSize", 14);
xlim([0, 5]);
ylim([-6000, 4100]);
saveas(gcf, "images/fourth_part_out3_u3.png");
%%

figure;
plot(t3, e_3(1,:), "LineWidth", 3);
hold on;
plot(t3, e_3(2, :), "LineWidth", 3);
plot(t3, e_3(3, :), "LineWidth", 3);
plot(t3, e_3(4, :), "LineWidth", 3);
title("Ошибка оценки e_3(t) при K_3 и L_3", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("e_3", "FontSize", 15);
legend("e_3(1)", "e_3(2)", "e_3(3)", "e_3(4)", "FontSize", 13, "Location", "northeast");
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/fourth_part_out3_e_3.png");

figure;
plot(t3, x_real_3(:, 1), "LineWidth", 3);
hold on;
plot(t3, x_est_3(:, 1), "LineWidth", 3, "LineStyle", "-.");
title("x_1(t) при K_3 и L_3", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_1", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "northeast");
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/fourth_part_out3_x_1.png");

figure;
plot(t3, x_real_3(:, 2), "LineWidth", 3);
hold on;
plot(t3, x_est_3(:, 2), "LineWidth", 3, "LineStyle", "-.");
title("x_2(t) при K_3 и L_3", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_2", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "northeast");
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/fourth_part_out3_x_2.png");

figure;
plot(t3, x_real_3(:, 3), "LineWidth", 3);
hold on;
plot(t3, x_est_3(:, 3), "LineWidth", 3, "LineStyle", "-.");
title("x_3(t) при K_3 и L_3", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_3", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "northeast");
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/fourth_part_out3_x_3.png");

figure;
plot(t3, x_real_3(:, 4), "LineWidth", 3);
hold on;
plot(t3, x_est_3(:, 4), "LineWidth", 3, "LineStyle", "-.");
title("x_4(t) при K_3 и L_3", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_4", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "northeast");
set(gca, "FontSize", 14);
xlim([0, 5]);
saveas(gcf, "images/fourth_part_out3_x_4.png");



figure;
plot(t4, u4, "LineWidth", 3);
title("Управление u_4(t) при K_4 и L_4", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("u_4", "FontSize", 15);
set(gca, "FontSize", 14);
xlim([0, 3]);
saveas(gcf, "images/fourth_part_out4_u4.png");

figure;
plot(t4, e_4(1,:), "LineWidth", 3);
hold on;
plot(t4, e_4(2, :), "LineWidth", 3);
plot(t4, e_4(3, :), "LineWidth", 3);
plot(t4, e_4(4, :), "LineWidth", 3);
title("Ошибка оценки e_4(t) при K_4 и L_4", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("e_4", "FontSize", 15);
legend("e_4(1)", "e_4(2)", "e_4(3)", "e_4(4)", "FontSize", 13, "Location", "southeast");
set(gca, "FontSize", 14);
xlim([0, 3]);
saveas(gcf, "images/fourth_part_out4_e_4.png");

figure;
plot(t4, x_real_4(:, 1), "LineWidth", 3);
hold on;
plot(t4, x_est_4(:, 1), "LineWidth", 3, "LineStyle", "-.");
title("x_1(t) при K_4 и L_4", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_1", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "northeast");
set(gca, "FontSize", 14);
xlim([0, 3]);
saveas(gcf, "images/fourth_part_out4_x_1.png");

%%
figure;
plot(t4, x_real_4(:, 2), "LineWidth", 3);
hold on;
plot(t4, x_est_4(:, 2), "LineWidth", 3, "LineStyle", "-.");
title("x_2(t) при K_4 и L_4", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_2", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "northeast");
set(gca, "FontSize", 14);
xlim([0, 3]);
ylim([-1.6, 2]);
saveas(gcf, "images/fourth_part_out4_x_2.png");
%%

figure;
plot(t4, x_real_4(:, 3), "LineWidth", 3);
hold on;
plot(t4, x_est_4(:, 3), "LineWidth", 3, "LineStyle", "-.");
title("x_3(t) при K_4 и L_4", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_3", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "northeast");
set(gca, "FontSize", 14);
xlim([0, 3]);
saveas(gcf, "images/fourth_part_out4_x_3.png");

figure;
plot(t4, x_real_4(:, 4), "LineWidth", 3);
hold on;
plot(t4, x_est_4(:, 4), "LineWidth", 3, "LineStyle", "-.");
title("x_4(t) при K_4 и L_4", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("x_4", "FontSize", 15);
legend("Истинное", "Оценка", "FontSize", 13, "Location", "northeast");
set(gca, "FontSize", 14);
xlim([0, 3]);
saveas(gcf, "images/fourth_part_out4_x_4.png");

close all;

%%
u1max = max(abs(u1))
u2max = max(abs(u2))
u3max = max(abs(u3))
u4max = max(abs(u4))

xa1max = max(abs(x_real_1(:, 1)))
xa2max = max(abs(x_real_2(:, 1)))
xa3max = max(abs(x_real_3(:, 1)))
xa4max = max(abs(x_real_4(:, 1)))

xp1max = max(abs(x_real_1(:, 3)))
xp2max = max(abs(x_real_2(:, 3)))
xp3max = max(abs(x_real_3(:, 3)))
xp4max = max(abs(x_real_4(:, 3)))