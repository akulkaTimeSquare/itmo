% Параметры системы
T = 0.2;

% Пять случаев для z_1 и z_2
cases = {
    [0.1, 0.7],           % Случай 1
    [-1.2, -0.4],         % Случай 2
    [0.1, 0.5],           % Случай 3
    [1.2j, -1.2j],        % Случай 4: z_{12} = ± 1.2 j
    [-0.8+0.7j, -0.8-0.7j] % Случай 5: z_{12} = -0.8 ± 0.7j
};

fprintf('Решение системы уравнений для различных z_1 и z_2:\n');
fprintf('Система:\n');
fprintf('  z_1 + z_2 = 2 - k_1*T^2/2 - k_2*T\n');
fprintf('  z_1*z_2 = 1 - k_1*T^2/2 - k_2*T + k_1*T^3/2\n\n');

% Решаем для каждого случая
for i = 1:length(cases)
    z1 = cases{i}(1);
    z2 = cases{i}(2);
    
    % Вычисляем сумму и произведение корней
    s = z1 + z2;  % сумма
    p = z1 * z2;  % произведение
    
    % Решаем систему аналитически:
    % Из уравнений получаем:
    % k_1 = 2*(1 - s + p) / T^3
    % k_2 = (2 - s - k_1*T^2/2) / T
    
    k1 = (1 - s + p) * 2 / T^3;
    k2 = (2-s) / T - (1 - s + p)/T^2;
    
    % Формируем вектор K
    K = [k1; k2];
    
    % Вывод результатов
    fprintf('Случай %d:\n', i);
    if isreal(z1) && isreal(z2)
        fprintf('  z_1 = %.1f, z_2 = %.1f\n', z1, z2);
    else
        fprintf('  z_1 = %s, z_2 = %s\n', num2str(z1), num2str(z2));
    end
    
    if isreal(K)
        fprintf('  K_%d = [%.4f; %.4f]\n', i, K(1), K(2));
    else
        fprintf('  K_%d = [%s; %s]\n', i, num2str(K(1)), num2str(K(2)));
    end
    
    % Проверка: подставляем обратно в систему
    check1 = 2 - k1*T^2/2 - k2*T;
    check2 = 1 - k1*T^2/2 - k2*T + k1*T^3/2;
    
    fprintf('  Проверка:\n');
    fprintf('    z_1 + z_2 = %s (должно быть %s)\n', num2str(s), num2str(check1));
    fprintf('    z_1*z_2 = %s (должно быть %s)\n', num2str(p), num2str(check2));
    fprintf('    Погрешность: %.2e, %.2e\n\n', abs(s - check1), abs(p - check2));
end

