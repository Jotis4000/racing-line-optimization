function laptime = calcLapTime(line)

    delta_pos = diff(line.geom);
    dx = hypot(delta_pos(:,1), delta_pos(:,2));
    dt = dx ./ line.V(1:(length(dx)));
    laptime = sum(dt);

end