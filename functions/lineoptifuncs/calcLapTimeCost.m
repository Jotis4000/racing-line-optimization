function cost = calcLapTimeCost(alpha_ctrl, s_ctrl, s_full, track, par)
    
    addpath("functions\")

    N = length(s_full);

    % Calculate line parameters the same as in the geometriic optimization
    alpha_full = makima(s_ctrl, alpha_ctrl, s_full);
    nx = track.vecleft(:,1) ./ track.vecmag;
    ny = track.vecleft(:,2) ./ track.vecmag;
    X_race = track.m(:,1) + alpha_full .* nx;
    Y_race = track.m(:,2) + alpha_full .* ny;
    
    dx = gradient(X_race);
    dy = gradient(Y_race);
    ddx = gradient(dx);
    ddy = gradient(dy);
    
    % Calculate curvature
    kappa = (dx .* ddy - dy .* ddx) ./ ((dx.^2 + dy.^2).^(3/2));
    % Calculate length
    ds = sqrt(dx.^2 + dy.^2);

    % Calculate the maximum cornering speed at all 900 points
    % v_profile = sqrt((track.mu * 9.81) ./ abs(kappa)); 
    % 
    % % Cap the maximum speed to your engine's top speed
    % v_profile = min(v_profile, par.Vmax);
    % 
    % for i = (N-1):-1:1
    %     % The max speed we can be at point i, to safely brake down to point i+1
    %     v_brake = sqrt(v_profile(i+1)^2 + 2 * par.decelmax*9.81 * ds(i));
    % 
    %     % The profile becomes whichever is slower: the grip limit, or the braking limit
    %     v_profile(i) = min(v_profile(i), v_brake);
    % end
    % 
    % for i = 2:N
    %     % The max speed we can reach at point i, accelerating from point i-1
    %     v_accel = sqrt(v_profile(i-1)^2 + 2 * par.accelmax*9.81 * ds(i-1));
    % 
    %     % The profile becomes whichever is slower
    %     v_profile(i) = min(v_profile(i), v_accel);
    % end

    v2_profile = (track.mu * 9.81) ./ abs(kappa); % Maximum speed per curvature of track
    v2_profile = min(v2_profile, par.Vmax^2); % At straight, max possible speed is Vmax
    
    % Deltas assuming constant maximum accel. and braking
    delta_v2_brake = 2 * par.decelmax*9.81 .* ds;
    delta_v2_accel = 2 * par.accelmax*9.81 .* ds;
    
    % 3. Backward Integration
    for i = (N-1):-1:1
        v2_profile(i) = min(v2_profile(i), v2_profile(i+1) + delta_v2_brake(i));
    end
    
    % Forward Integration
    for i = 2:N
        v2_profile(i) = min(v2_profile(i), v2_profile(i-1) + delta_v2_accel(i-1));
    end
    
    v_profile = sqrt(v2_profile);

    cost = sum(ds./v_profile);

end

