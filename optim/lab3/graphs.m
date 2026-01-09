% Constant
s6 = sqrt(6);

% A_i coefficients
A1 = 2*exp(s6+3) + 2*exp(s6+9) + 2*exp(s6+15);
A2 = -6*exp(3-s6) + 6*exp(9-3*s6) - 6*exp(15-s6);
A3 = s6 * ( ...
    -2*exp(3-s6) ...
    -2*exp(9-s6) ...
    +3*exp(9-3*s6) ...
    -2*exp(15-s6) ...
    +3*exp(s6+9) );
A = A1 + A2 + A3;

% B_i coefficients
B1 = 18*exp(s6+3) - 30*exp(3-s6) - 12*exp(9-s6);
B2 = 36*exp(9-3*s6) + 6*exp(15-s6) - 18*exp(15-3*s6);
B3 = ...
    -12*s6*exp(3-s6) ...
    -6*s6*exp(9-s6) ...
    +15*s6*exp(9-3*s6) ...
    -6*s6*exp(15-3*s6) ...
    +6*s6*exp(s6+3) ...
    +3*s6*exp(s6+9);
B = B1 + B2 + B3;

% Denominator D
D1 = s6*exp(6) ...
    -6*exp(6-2*s6) ...
    -3*s6*exp(6-2*s6) ...
    +2*s6 + 6;
D2 = 2*exp(6) + 2*exp(12) ...
    -3*exp(6-2*s6) ...
    -3*exp(2*s6+6) + 2;
D = D1 * D2;

% Control law u(t)
e = -0.0001;
% J = 6829.24;
u = @(t) ...
    (1 + e)*(180*exp(-t*(s6-3))/D)*A ...
  - (1 + e)*(60*exp(t*(s6+3))/D)*B;

% System dynamics
f = @(t,x) [ ...
    x(2);
   -3*x(1) - 6*x(2) + u(t) ];

% Simulation
tspan = [0 1];
x0 = [0; 0];

opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,x] = ode45(f, tspan, x0, opts);

x_final = x(end,:).';

% Output
disp('Final state x(1):');
disp(x_final);

disp('Target state:');
disp([10; 0]);

disp('Error:');
disp(x_final - [10; 0]);

%% Plots x
figure;
plot(t, x(:,1), 'LineWidth', 2); hold on;
plot(t, x(:,2), 'LineWidth', 2)
hold off;
title('Переменные состояния объекта');
xlabel('t');
ylabel('x(t)');
grid on;
legend("x_1(t)", "x_2(t)");
saveas(gcf, "images/x1.png");

%% Plots u
figure;
tt = linspace(0, 1, 1000);
plot(tt, arrayfun(u,tt), 'LineWidth', 2)
title('Оптимальное управление');
xlabel('t');
ylabel('u(t)');
grid on;
saveas(gcf, "images/u1.png");

%% Performance index J = ∫ u^2 dt
J = integral(@(t) u(t).^2, 0, 1, ...
             'RelTol',1e-10,'AbsTol',1e-12);

disp('Performance index J = ∫ u^2 dt:');
disp(J);

%% J
t_fine = linspace(0, 1, 50000);
u_fine = arrayfun(u, t_fine);
u_sq_fine = u_fine.^2;

J_cum = cumtrapz(t_fine, u_sq_fine);
J_end_num = J_cum(end);
fprintf('Контроль накопления: J(1) = %.6f, J (из integral) = %.6f, разница = %.2e\n', ...
        J_end_num, J, abs(J_end_num - J));

% --- График накопленного критерия ---
figure;
plot(t_fine, J_cum, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]); % оранжевый
title('Накопленный критерий качества');
xlabel('t');
ylabel('J(t)');
grid on;
yticks([0 1000 2000 3000 4000 5000 6000 6829]);
saveas(gcf, "images/J_cum1.png");