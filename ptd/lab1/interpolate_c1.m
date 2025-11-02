function pts = interpolate_c1(points, r, step)
    % C1 непрерывность через сглаживание дугами окружностей в углах
    pts = [];
    p_exit_last = points(1, :);

    for i = 2:(size(points, 1) - 1)
        p1 = points(i - 1, :);
        p2 = points(i, :);
        p3 = points(i + 1, :);

        corner_angle = compute_angle(p1, p2, p3);

        if abs(corner_angle - pi) < 1e-6
            pts = [pts; interpolate_c0([p_exit_last; p2], step)];
            p_exit_last = p2;
            continue;
        end

        dir1 = p2 - p1; dir1 = dir1 / norm(dir1);
        dir2 = p3 - p2; dir2 = dir2 / norm(dir2);

        sign_dir = sign(det([dir1; dir2]));  % аналог np.cross для 2D
        half_corner = corner_angle / 2;

        d = r / abs(tan(half_corner));
        dc = r / abs(sin(half_corner));
        phi1 = atan2(dir1(2), dir1(1));

        cc = rot_mat(pi - sign_dir * half_corner + phi1) * [dc; 0] + p2(:);

        p_entry = p2 - dir1 * d;
        p_exit  = p2 + dir2 * d;

        n_points = max(8, floor(corner_angle / step * r / 2));

        arc = interpolate_arc(cc', r, p_entry, p_exit, sign_dir, n_points);

        pts = [pts; interpolate_c0([p_exit_last; p_entry], step)];
        pts = [pts; arc];

        p_exit_last = p_exit;
    end

    pts = [pts; interpolate_c0([p_exit_last; points(end, :)], step)];
end
