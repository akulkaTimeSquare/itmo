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

G1 = [-0.5  0 0 0; 0 -0.75 0 0; 0 0 -1 0; 0 0 0 -1.25];
G2 = [-0.5  0 0 0; 0 -0.75 0 0; 0 0 -1 0; 0 0 0 -1.25];
G3 = [-2 0 0 0; 0 -2.25 0 0; 0 0 -2.8 0; 0 0  0 -2.75];
G4 = [-2 0 0 0; 0 -2.25 0 0; 0 0 -2.8 0; 0 0  0 -2.75];

Gl1 = [-3.25 0 0 0; 0 -4 0 0; 0 0 -3.4 0; 0 0  0 -3.8];
Gl2 = [-8  0 0 0; 0 -7.5 0 0; 0 0 -8.25 0; 0 0 0 -8.5];
Gl3 = [-3.25 0 0 0; 0 -4 0 0; 0 0 -3.4 0; 0 0  0 -3.8];
Gl4  = [-8  0 0 0; 0 -7.5 0 0; 0 0 -8.25 0; 0 0 0 -8.5];

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

er1 = eig(A + B*K1)
er2 = eig(A + B*K2)
er3 = eig(A + B*K3)
er4 = eig(A + B*K4)

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

Yl = [1 1; 1 0; 1 1; 0 1];
V1 = [Yl Gl1*Yl Gl1^2*Yl Gl1^3*Yl];
V2 = [Yl Gl2*Yl Gl2^2*Yl Gl2^3*Yl];
V3 = [Yl Gl3*Yl Gl3^2*Yl Gl3^3*Yl];
V4 = [Yl Gl4*Yl Gl4^2*Yl Gl4^3*Yl];
display(rank(V1))
display(rank(V2))
display(rank(V3))
display(rank(V4))

Q1 = sylvester(Gl1, -A, Yl*C);
L1 = inv(Q1)*Yl;
L1T = transpose(L1);

Q2 = sylvester(Gl2, -A, Yl*C);
L2 = inv(Q2)*Yl;
L2T = transpose(L2);

Q3 = sylvester(Gl3, -A, Yl*C);
L3 = inv(Q3)*Yl;
L3T = transpose(L3);

Q4 = sylvester(Gl4, -A, Yl*C);
L4 = inv(Q4)*Yl;
L4T = transpose(L4);
disp("L1T:");
disp(L1T)
disp("L2T:");
disp(L2T)
disp("L3T:");
disp(L3T)
disp("L4T:");
disp(L4T)

el1 = eig(A + L1*C)
el2 = eig(A + L2*C)
el3 = eig(A + L3*C)
el4 = eig(A + L4*C)

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
xlim([0, 14]);
saveas(gcf, "images/third_part_reg_out_u1.png");

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
xlim([0, 5]);
saveas(gcf, "images/third_part_reg_out_e_1.png");

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
xlim([0, 14]);
saveas(gcf, "images/third_part_reg_out1_x_1.png");

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
xlim([0, 14]);
ylim([-0.7, 0.4]);
saveas(gcf, "images/third_part_reg_out1_x_2.png");

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
xlim([0, 14]);
saveas(gcf, "images/third_part_reg_out1_x_3.png");

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
xlim([0, 14]);
saveas(gcf, "images/third_part_reg_out1_x_4.png");



figure;
plot(t2, u2, "LineWidth", 3);
title("Управление u_2(t) при K_2 и L_2", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("u_2", "FontSize", 15);
set(gca, "FontSize", 14);
xlim([0, 14]);
saveas(gcf, "images/third_part_reg_out2_u2.png");

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
xlim([0, 3]);
ylim([-0.22, 0.1]);
saveas(gcf, "images/third_part_reg_out2_e_2.png");

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
xlim([0, 14]);
ylim([-0.32, 0.1]);
saveas(gcf, "images/third_part_reg_out2_x_1.png");

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
xlim([0, 14]);
ylim([-0.45, 0.2]);
saveas(gcf, "images/third_part_reg_out2_x_2.png");

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
xlim([0, 14]);
saveas(gcf, "images/third_part_reg_out2_x_3.png");

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
xlim([0, 14]);
saveas(gcf, "images/third_part_reg_out2_x_4.png");


figure;
plot(t3, u3, "LineWidth", 3);
title("Управление u_3(t) при K_3 и L_3", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("u_3", "FontSize", 15);
set(gca, "FontSize", 14);
xlim([0, 10]);
saveas(gcf, "images/third_part_reg_out3_u3.png");

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
saveas(gcf, "images/third_part_reg_out3_e_3.png");

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
xlim([0, 10]);
saveas(gcf, "images/third_part_reg_out3_x_1.png");

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
xlim([0, 10]);
ylim([-0.55, 0.35]);
saveas(gcf, "images/third_part_reg_out3_x_2.png");

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
xlim([0, 10]);
saveas(gcf, "images/third_part_reg_out3_x_3.png");

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
xlim([0, 10]);
saveas(gcf, "images/third_part_reg_out3_x_4.png");



figure;
plot(t4, u4, "LineWidth", 3);
title("Управление u_4(t) при K_4 и L_4", "FontSize", 18);
grid on;
xlabel("t", "FontSize", 15);
ylabel("u_4", "FontSize", 15);
set(gca, "FontSize", 14);
xlim([0, 6]);
saveas(gcf, "images/third_part_reg_out4_u4.png");

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
ylim([-0.22, 0.1]);
saveas(gcf, "images/third_part_reg_out4_e_4.png");

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
xlim([0, 6]);
saveas(gcf, "images/third_part_reg_out4_x_1.png");

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
xlim([0, 6]);
saveas(gcf, "images/third_part_reg_out4_x_2.png");

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
xlim([0, 6]);
saveas(gcf, "images/third_part_reg_out4_x_3.png");

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
xlim([0, 6]);
ylim([-0.45, 0.2]);
saveas(gcf, "images/third_part_reg_out4_x_4.png");

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