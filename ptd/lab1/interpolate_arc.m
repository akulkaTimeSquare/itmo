function pts = interpolate_arc(center, r, p_entry, p_exit, sign_dir, pts_num)
    start_angle = atan2(p_entry(2) - center(2), p_entry(1) - center(1));
    end_angle   = atan2(p_exit(2)  - center(2), p_exit(1)  - center(1));

    % Обеспечить правильное направление вращения
    if sign_dir > 0 && end_angle < start_angle
        end_angle = end_angle + 2*pi;
    elseif sign_dir < 0 && end_angle > start_angle
        end_angle = end_angle - 2*pi;
    end

    angles = linspace(start_angle, end_angle, pts_num);
    x = center(1) + r * cos(angles);
    y = center(2) + r * sin(angles);
    pts = [x(:), y(:)];
end
