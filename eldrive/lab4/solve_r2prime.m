function results = solve_r2prime_with_outputs(r1, lambda)
% Решение уравнения mu_m(r2') = lambda и вывод промежуточных параметров
%
% Пример запуска:
%   res = solve_r2prime_with_outputs(0.5, 2.6);

%% Данные двигателя
P_n   = 7500;        % Вт
n_n   = 1450;        % об/мин
eta   = 0.875;
cosphi = 0.85;
U_n   = 380;         % В (линейное)
I_n   = 15.3;        % A
M_n   = 49.4;        % Нм
k_s   = 2.1;
k_sI  = 7.0;
lambda_max = 2.8;    % относит. макс. момент
J     = 0.032;       % кгм^2
mass  = 70;          % кг
Fs    = 1.15;
m1    = 3;
f_s   = 50;          % Гц
z_p = 2;

U1N = U_n / sqrt(3);     % фазное напряжение
I1N = I_n;
omega1 = 2*pi*f_s;
n1 = 1500;               % синхронная скорость (об/мин, 4 полюса)
s_n = 1 - n_n/n1;

%% Локальная функция невязки
    function res = residual(r2p)
        a = r1 / r2p;
        D = lambda^2 - a;
        if D <= 0
            res = 1e6; return;
        end
        s_m = s_n * (lambda + sqrt(D));
        if s_m <= 0, res = 1e6; return; end
        term = (r2p/s_m)^2 - r1^2;
        if term <= 0, res = 1e6; return; end
        x_ks = sqrt(term);
        denom_b = (r1 + r2p/s_n)^2 + x_ks^2;
        b = x_ks/denom_b;
        sinphi = sqrt(1 - cosphi^2);
        denom_xm = (I1N*sinphi/U1N) - b;
        if denom_xm <= 0, res = 1e6; return; end
        x_m = 1/denom_xm;
        I2p = U1N / sqrt((r1 + r2p/s_m)^2 + x_ks^2);
        mu_m = (m1*z_p*abs(I2p)*r2p)/(omega1*s_m*M_n);
        res = mu_m - lambda;
    end

%% Решение уравнения
r_guess = max(1e-3, r1);
options = optimset('TolX',1e-9,'Display','off');
r2prime = fsolve(@(x) residual(x), r_guess, options);

%% Пересчёт параметров при найденном r2'
a = r1 / r2prime;
s_m = s_n * (lambda + sqrt(lambda^2 - a));
x_ks = sqrt((r2prime/s_m)^2 - r1^2);
b = x_ks / ((r1 + r2prime/s_n)^2 + x_ks^2);
x_s1 = x_ks/2;   % = x'_s2
sinphi = sqrt(1 - cosphi^2);
x_m = 1 / ((I1N*sinphi/U1N) - b);
I2p = U1N / sqrt((r1 + r2prime/s_m)^2 + x_ks^2);

%% Вывод
fprintf("r2'  = %.2f Ом\n", r2prime);
fprintf("a    = %.2f\n", a);
fprintf("s_m  = %.2f\n", s_m);
fprintf("b    = %.4f\n", b);
fprintf("xs1  = xs2 = %.2f\n", x_s1);
fprintf("xm   = %.2f\n", x_m);
fprintf("I2'  = %.2f A\n", I2p);

%% Сохранение в структуру
results = struct('r2prime',r2prime,'a',a,'s_m',s_m,'b',b,...
                 'xs1',x_s1,'xm',x_m,'I2p',I2p);
end

res = solve_r2prime_with_outputs(1.14, 2.8);