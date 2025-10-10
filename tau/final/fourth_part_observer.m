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

x0 = [0.05; 0.1; 0.1; -0.02];

ar = 1;

cvx_begin sdp quiet
variable P(4, 4) symmetric
variable Y(1, 4)
P > 0.0001*eye(4);
P*A' + A*P + 2*ar*P + Y'*B' + B*Y <= 0;
cvx_end

K = Y*inv(P);
er = eig(A + B*K)

al = 1.5;

cvx_begin sdp quiet
variable Q(4, 4) symmetric
variable Y(4, 2)
Q > 0.0001*eye(4);
A'*Q + Q*A + 2*al*Q + C'*Y' + Y*C <= 0;
cvx_end

L = Q \ Y;
el = eig(A + L*C)

function printMatrix(Min, prec)
    fmt = ['%0.', num2str(prec), 'f  '];   % формат вывода
    for i = 1:size(Min,1)
        fprintf(fmt, Min(i,:));
        fprintf('\n');
    end
end

disp("K:")
printMatrix(K, 4)
disp("L:")
printMatrix(L', 4)


% Nonlinear system simulation with observer'
dt = 0.0001;
t_values = 0:dt:15;

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
x_hat0 = [0; 0; 0; 0];
z0 = [x0; x_hat0];

% Solve for different observers
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

[t, sol] = ode45(@(t,z) nonlinear_system_observer(t, z, L, A, B, C, K, M, m, l, g), t_values, z0, options);

% Extract real and estimated states
x_real = sol(:, 1:4);
x_est = sol(:, 5:8);

% Calculate errors
e = x_real - x_est;

figure;
plot(t, e(:,1), "LineWidth", 2);
hold on;
plot(t, e(:,2), "LineWidth", 2);
plot(t, e(:,3), "LineWidth", 2);
plot(t, e(:,4), "LineWidth", 2);
title("Ошибка оценки e при \alpha_L = 1.5");
grid on;
xlabel("t");
ylabel("e");
xlim([0, 5]);
legend("e(1)", "e(2)", "e(3)", "e(4)", "Location", "northeast");
saveas(gcf, "images/fourth_part_observer_e.png");

figure;
plot(t, x_real(:,1), "LineWidth", 2, "Clipping", "on"); hold on;
plot(t, x_est(:,1), "LineWidth", 2, LineStyle="-.");
title("x_1(t) при \alpha_L = 1.5");
grid on;
xlabel("t");
ylabel("x_1");
xlim([0, 5]);
legend("Истинное", "Оценка", "Location", "northeast");
saveas(gcf, "images/fourth_part_observer_x_1.png");

figure;
plot(t, x_real(:,2), "LineWidth", 2, "Clipping", "on"); hold on;
plot(t, x_est(:,2), "LineWidth", 2, LineStyle="-.");
title("x_2(t) при \alpha_L = 1.5");
grid on;
xlabel("t");
ylabel("x_2");
xlim([0, 5]);
ylim([-0.5 0.25]);
legend("Истинное", "Оценка", "Location", "northeast");
saveas(gcf, "images/fourth_part_observer_x_2.png");

figure;
plot(t, x_real(:,3), "LineWidth", 2, "Clipping", "on"); hold on;
plot(t, x_est(:,3), "LineWidth", 2, LineStyle="-.");
title("x_3(t) при \alpha_L = 1.5");
grid on;
xlabel("t");
ylabel("x_3");
xlim([0, 5]);
legend("Истинное", "Оценка", "Location", "northeast");
saveas(gcf, "images/fourth_part_observer_x_3.png");

figure;
plot(t, x_real(:,4), "LineWidth", 2, "Clipping", "on"); hold on;
plot(t, x_est(:,4), "LineWidth", 2, LineStyle="-.");
title("x_4(t) при \alpha_L = 1.5");
grid on;
xlabel("t");
ylabel("x_4");
xlim([0, 5]);
legend("Истинное", "Оценка", "Location", "northeast");
saveas(gcf, "images/fourth_part_observer_x_4.png");