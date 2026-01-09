function dydt = optimal_ode_system(t, y)
% Описание полной системы ОДУ для BVP решателя.
% Вектор y = [x1; x2; f1; f2]

% Управление u = f2/2
u = y(4) / 2;

dydt = [
    y(2);                    % dx1/dt = x2
    -3*y(1) - 6*y(2) + u;    % dx2/dt = -3x1 - 6x2 + u
    3*y(4);                  % df1/dt = 3*f2
    -y(3) + 6*y(4)           % df2/dt = -f1 + 6*f2
];
end
