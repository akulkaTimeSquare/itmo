function angle = compute_angle(p1, p2, p3)
    v1 = p1 - p2;
    v2 = p3 - p2;
    a = norm(v1);
    b = norm(v2);

    if a == 0 || b == 0
        angle = 0; return;
    end

    cos_angle = max(min(dot(v1, v2) / (a * b), 1.0), -1.0);
    angle = acos(cos_angle);

    if angle < 1e-6
        angle = 0; return;
    end

    if abs(angle - pi) < 1e-5
        angle = pi; return;
    end
end
