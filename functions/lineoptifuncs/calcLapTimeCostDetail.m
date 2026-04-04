function cost = calcLapTimeCostDetail(alpha_ctrl, s_ctrl, s_full, track, par)
    
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

    g = 9.81;
    
    Lconst = 0.5 * par.rho * par.CLA; % Downforce constant
    Dconst = 0.5 * par.rho * par.CDA; % Drag constant
    
    % par.Vmax = sqrt(par.F_engine/(par.CDA*0.5*par.rho));

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
        % 1. What is the speed at the NEXT point?
        v2_next = v2_profile(i+1);
        
        % 2. Calculate Total Normal Force at that speed
        F_normal = par.m * g + (Lconst * v2_next);
        
        % 3. Calculate Max Braking Force (Grip limit)
        % Note: If you add an aerodynamic DRAG term (CDA), you add it here!
        F_brake_total = (track.mu * F_normal) + (Dconst * v2_next);
        % F_brake_total = track.mu * F_normal; 
        
        % 4. Convert to deceleration and update
        a_brake = F_brake_total / par.m;
        delta_v2 = 2 * a_brake * ds(i);
        
        v2_profile(i) = min(v2_profile(i), v2_next + delta_v2);
    end
    
    %% PASS 3: Forward Integration (Acceleration Zones)
    for i = 2:N
        v2_prev = v2_profile(i-1);
        
        % We need true velocity to calculate engine thrust
        v_prev = sqrt(v2_prev); 
        
        % 1. Calculate Aerodynamic Forces
        F_downforce = Lconst * v2_prev;
        F_drag      = Dconst * v2_prev;
        
        % 2. Calculate Total Normal Force & Available Grip
        F_normal = par.m * 9.81 + F_downforce;
        F_grip_avail = track.mu * F_normal;
        
        % 3. Calculate Engine Thrust (The Power Model)
        if v_prev < 1.0 % Prevent divide-by-zero at standstill
            F_engine = par.F_engine_max;
        else
            F_engine = min(par.P_engine / v_prev, par.F_engine_max);
        end
        
        % 4. Calculate Net Engine Thrust
        F_thrust_net = F_engine - F_drag;
        
        % 5. Actual acceleration force (Traction Limited vs Power Limited)
        F_accel = min(F_grip_avail, F_thrust_net);
        
        % 6. Convert to acceleration and update
        a_accel = F_accel / par.m;
        delta_v2 = 2 * a_accel * ds(i-1);
        
        v2_profile(i) = min(v2_profile(i), v2_prev + delta_v2);
    end
    
    v_profile = sqrt(v2_profile);
    cost = 10*sum(ds./v_profile);

end

