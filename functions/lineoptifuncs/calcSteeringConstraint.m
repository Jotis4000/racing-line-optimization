function [c, ceq] = calcSteeringConstraint(alpha_ctrl, s_ctrl, s_full, track, par)
    % 1. Rebuild the full racing line from control points
    alpha_full = makima(s_ctrl, alpha_ctrl, s_full);
    nx = track.vecleft(:,1) ./ track.vecmag;
    ny = track.vecleft(:,2) ./ track.vecmag;
    X_race = track.m(:,1) + alpha_full .* nx;
    Y_race = track.m(:,2) + alpha_full .* ny;
    
    % 2. Calculate derivatives to find curvature
    dx = gradient(X_race);
    dy = gradient(Y_race);
    ddx = gradient(dx);
    ddy = gradient(dy);
    
    kappa = (dx .* ddy - dy .* ddx) ./ ((dx.^2 + dy.^2).^(3/2));
    
    % 3. Define the maximum allowable curvature
    kappa_max = 1 / par.R_min;
    
    % 4. Inequality constraint: c <= 0
    % We want abs(kappa) - kappa_max <= 0 at EVERY point
    c = abs(kappa) - kappa_max; 
    
    % 5. No equality constraints needed for this
    ceq = []; 
end