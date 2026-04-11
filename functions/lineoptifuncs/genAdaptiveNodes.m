function s_ctrl = genAdaptiveNodes(s_full, track_m, n_var, corner_weight, smoothing_window)
    % GENADAPTIVENODES Automatically clusters control points around corners.
    
    % 1. Calculate raw curvature (using standard central differences)
    dx = gradient(track_m(:,1));
    dy = gradient(track_m(:,2));
    ddx = gradient(dx);
    ddy = gradient(dy);
    kappa_raw = (dx .* ddy - dy .* ddx) ./ ((dx.^2 + dy.^2).^(3/2));

    % 2. Smooth the curvature to prevent stuttering on micro-apexes
    kappa_smooth = smoothdata(abs(kappa_raw), 'gaussian', smoothing_window); 

    % 3. Normalize and square to aggressively target the tightest hairpins
    kappa_norm = kappa_smooth / max(kappa_smooth);
    kappa_focused = kappa_norm .^ 2; 

    % 4. Define node density and integrate
    node_density = 1 + (corner_weight .* kappa_focused);
    cum_density = cumtrapz(s_full, node_density);

    % 5. Sample evenly in the density domain 
    cum_density_samples = linspace(0, cum_density(end), n_var)';
    
    % 6. Map back to the physical track distance
    s_ctrl = interp1(cum_density, s_full, cum_density_samples, 'pchip');
end