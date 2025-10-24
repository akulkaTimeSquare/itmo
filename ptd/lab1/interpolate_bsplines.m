function pts = interpolate_bsplines(waypoints, params)
    n = size(waypoints, 1);
    if n < 2
        pts = waypoints;
        return;
    end

    k = params.bspline_degree;      % степень
    order = k + 1;                  % порядок

    if n <= order
        order = n - 1;
        k = order - 1;
    end

    % --- Узловой вектор (clamped, с повтором концов) ---
    knots = [zeros(1, order), 1:1:n - k, (n-k)*ones(1, order)];

    % --- Параметризация по t ---
    t = 0:params.resolution:n-k;
    pts = zeros(length(t), size(waypoints, 2));

    for i = 1:n
        Ni = bspline_basis_vectorized(t, i, k, knots)';
        pts = pts + Ni .* waypoints(i, :);
    end

    % Гарантируем попадание в последнюю точку
    pts(end, :) = waypoints(end, :);
end


function N = bspline_basis(t, i, k, knots)
    % B-spline базис (рекурсивное определение Кокса — де Бура)
    if k == 0
        N = double(knots(i) <= t & t < knots(i+1));
        return;
    end

    leftDen = knots(i+k) - knots(i);
    rightDen = knots(i+k+1) - knots(i+1);

    left = 0; right = 0;
    if leftDen ~= 0
        left = ((t - knots(i)) ./ leftDen) .* bspline_basis(t, i, k-1, knots);
    end
    if rightDen ~= 0
        right = ((knots(i+k+1) - t) ./ rightDen) .* bspline_basis(t, i+1, k-1, knots);
    end

    N = left + right;
end


function vals = bspline_basis_vectorized(ts, i, k, knots)
    vals = arrayfun(@(t) bspline_basis(t, i, k, knots), ts);
end