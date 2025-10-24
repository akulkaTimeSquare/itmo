function pts = interpolate_c2(waypoints, r, k, step)
    % C2 continuity: соединение прямых, кубических парабол и дуг
    pts = [];
    p_exit_last = waypoints(1, :);

    for i = 2:(size(waypoints, 1) - 1)
        a = waypoints(i - 1, :);
        b = waypoints(i, :);
        c = waypoints(i + 1, :);

        corner_angle = compute_angle(a, b, c);

        if abs(corner_angle - pi) < 1e-6
            pts = [pts; interpolate_c0([p_exit_last; b], step)];
            p_exit_last = b;
            continue;
        end

        half_corner = corner_angle / 2;
        d = r / abs(tan(half_corner));

        dir1 = b - a; dir1 = dir1 / norm(dir1);
        dir2 = c - b; dir2 = dir2 / norm(dir2);
        sign_dir = sign(det([dir1; dir2]));

        phi1 = atan2(dir1(2), dir1(1));
        phi2 = atan2(dir2(2), dir2(1));

        parab1_entry = b - dir1 * d;
        parab1_exit = rot_mat(phi1) * [1 / (6 * k * r); sign_dir * k / (6 * k * r)^3] + parab1_entry(:);
        parab1_exit = parab1_exit(:)';

        parab2_exit = b + dir2 * d;
        parab2_entry = rot_mat(phi2 - pi) * [1 / (6 * k * r); -sign_dir * k / (6 * k * r)^3] + parab2_exit(:);
        parab2_entry = parab2_entry(:)';

        h = norm(parab2_entry - parab1_exit);
        dc1 = sqrt(r^2 - (h / 2)^2);
        dc2 = sqrt(norm(b - parab2_entry)^2 - (h / 2)^2);
        dc = dc1 + dc2;

        cc = rot_mat(pi - sign_dir * half_corner + phi1) * [dc; 0] + b(:);
        cc = cc(:)';

        % Составление сегментов
        line = interpolate_c0([p_exit_last; parab1_entry], step);

        parab1 = interpolate_cubic_parabola(parab1_exit - parab1_entry, sign_dir * k, step, phi1);
        parab1 = parab1 + parab1_entry;

        parab2 = interpolate_cubic_parabola(parab2_exit - parab2_entry, -sign_dir * k, step, phi2 - pi);
        parab2 = flipud(parab2) + parab2_exit;

        n_points = max(20, floor(abs(corner_angle) / step * r / 2));
        arc = interpolate_arc(cc, r, parab1_exit, parab2_entry, sign_dir, n_points);

        pts = [pts; line; parab1; arc; parab2];
        p_exit_last = parab2_exit;
    end

    pts = [pts; interpolate_c0([p_exit_last; waypoints(end, :)], step)];
end

function pts = interpolate_cubic_parabola(p_exit, k, resolution, rotation)
    if nargin < 4
        rotation = 0;
    end

    R = rot_mat(rotation);
    tmp = R' * p_exit(:);
    xL_end = abs(tmp(1));
    xL = 0:resolution:xL_end;

    yL = k * (xL.^3);

    pts_local = [xL; yL];
    pts = (R * pts_local)';
end
