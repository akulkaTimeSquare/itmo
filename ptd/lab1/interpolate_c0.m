function pts = interpolate_c0(waypoints, step)
    % C0 continuity (piecewise linear)
    if size(waypoints, 1) < 2
        pts = waypoints;
        return;
    end

    pts = [];
    for i = 1:(size(waypoints, 1) - 1)
        p1 = waypoints(i, :);
        p2 = waypoints(i + 1, :);
        n = max(2, floor(norm(p2 - p1) / step)); % шаг по длине
        segment = [linspace(p1(1), p2(1), n)', linspace(p1(2), p2(2), n)'];
        pts = [pts; segment];
    end
end
