function cost = calcLapTimeCostDetailAdapt(alpha_ctrl, s_ctrl, s_full, track, par, splineType, smooth_weight)
    % calcLapTimeCostDetailAdapt - Calculates lap time cost with dynamic interpolation
    
    % --- DEFAULT ARGUMENTS ---
    if nargin < 6
        splineType = 'spline'; % Default interpolation
    end
    if nargin < 7
        smooth_weight = 0.0;   % Default to no smoothness penalty
    end

    N = length(s_full);

    % --- INTERPOLATION BLOCK ---
    switch lower(splineType)
        case 'makima'
            alpha_full = makima(s_ctrl, alpha_ctrl, s_full);
        case 'bspline'
            % Uses MATLAB Curve Fitting Toolbox B-Spline
            b_curve = spapi(4, s_ctrl, alpha_ctrl);
            alpha_full = fnval(b_curve, s_full);
        otherwise
            % Default to standard cubic spline
            alpha_full = spline(s_ctrl, alpha_ctrl, s_full);
    end

    % CRITICAL FIX: Ensure alpha_full is always a column vector 
    % so it doesn't crash when multiplied by nx and ny
    alpha_full = alpha_full(:);

    % Calculate line parameters
    nx = track.vecleft(:,1) ./ track.vecmag;
    ny = track.vecleft(:,2) ./ track.vecmag;
    X_race = track.m(:,1) + alpha_full .* nx;
    Y_race = track.m(:,2) + alpha_full .* ny;
    
    dx = gradient(X_race);
    dy = gradient(Y_race);
    ddx = gradient(dx);
    ddy = gradient(dy);
    
    % Calculate curvature and length
    kappa = (dx .* ddy - dy .* ddx) ./ ((dx.^2 + dy.^2).^(3/2));
    ds = sqrt(dx.^2 + dy.^2);

    g = 9.81;
    Lconst = 0.5 * par.rho * par.CLA; % Downforce constant
    Dconst = 0.5 * par.rho * par.CDA; % Drag constant

    %% PASS 1: The Aero Grip Ceiling
    denominator = (par.m .* abs(kappa)) - (track.mu * Lconst);
    v2_profile = zeros(N, 1);
    
    for i = 1:N
        if denominator(i) <= 0
            % Downforce outscales cornering demand. Flat out!
            v2_profile(i) = par.Vmax^2;
        else
            % Calculate limit, cap at vehicle top speed
            v2_limit = (track.mu * par.m * g) / denominator(i);
            v2_profile(i) = min(v2_limit, par.Vmax^2);
        end
    end
    
    %% PASS 2: Backward Integration (Braking Zones)
    for i = (N-1):-1:1
        v2_next = v2_profile(i+1);
        
        F_normal = par.m * g + (Lconst * v2_next);
        F_grip_total = track.mu * F_normal;
        F_lat = par.m * v2_next * abs(kappa(i+1));
        
        % max(..., 0) prevents imaginary numbers if slight numerical errors occur
        F_brake_avail = sqrt(max(F_grip_total^2 - F_lat^2, 0));
        F_brake_total = F_brake_avail + (Dconst * v2_next);
        
        a_brake = F_brake_total / par.m;
        delta_v2 = 2 * a_brake * ds(i);
        
        v2_profile(i) = min(v2_profile(i), v2_next + delta_v2);
    end
    
    %% PASS 3: Forward Integration (Acceleration Zones)
    for i = 2:N
        v2_prev = v2_profile(i-1);
        v_prev = sqrt(v2_prev); 
        
        F_downforce = Lconst * v2_prev;
        F_drag      = Dconst * v2_prev;
        
        F_normal = par.m * g + F_downforce;
        F_grip_total = track.mu * F_normal;
        
        F_lat = par.m * v2_prev * abs(kappa(i-1));
        F_grip_long_avail = sqrt(max(F_grip_total^2 - F_lat^2, 0));
        
        if v_prev < 1.0 
            F_engine = par.F_engine_max;
        else
            F_engine = min(par.P_engine / v_prev, par.F_engine_max);
        end
        F_thrust_net = F_engine - F_drag;
        
        F_accel = min(F_grip_long_avail, F_thrust_net);
        
        a_accel = F_accel / par.m;
        delta_v2 = 2 * a_accel * ds(i-1);
        
        v2_profile(i) = min(v2_profile(i), v2_prev + delta_v2);
    end
    
    v_profile = sqrt(v2_profile);
    
    % --- FINAL COST CALCULATION ---
    laptime_cost = 10 * sum(ds./v_profile);
    
    % Add the smoothness penalty based on the control points (wiggle penalty)
    alpha_wiggle = diff(alpha_ctrl, 2);
    smoothness_penalty = sum(alpha_wiggle.^2);
    
    cost = laptime_cost + (smooth_weight * smoothness_penalty);

end