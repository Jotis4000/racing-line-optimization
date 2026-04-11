% function cost = calcCurvatureCost(alpha_ctrl, s_ctrl, s_full, track)
% 
%     alpha_full = makima(s_ctrl, alpha_ctrl, s_full); % Spline is sampled back to full detail
% 
%     % 2. Get normalized normal vectors 
%     % (pointing left based on your genTrack function)
%     nx = track.vecleft(:,1) ./ track.vecmag;
%     ny = track.vecleft(:,2) ./ track.vecmag;
% 
%     % 3. Reconstruct the X, Y racing line
%     X_race = track.m(:,1) + alpha_full .* nx;
%     Y_race = track.m(:,2) + alpha_full .* ny;
% 
%     % 4. Calculate Curvature (kappa)
%     % First derivatives (velocity vectors)
%     dx = gradient(X_race);
%     dy = gradient(Y_race);
% 
%     % Second derivatives (acceleration vectors)
%     ddx = gradient(dx);
%     ddy = gradient(dy);
% 
%     % Curvature formula
%     kappa = (dx .* ddy - dy .* ddx) ./ ((dx.^2 + dy.^2).^(3/2));
% 
%     % 5. Objective: Minimize Sum of Squared Curvature
%     cost = 1000*sum(kappa.^2);
% 
% end

function cost = calcCurvatureCost(alpha_ctrl, s_ctrl, s_full, track, weight_length)
    
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
 
    ds_race = sqrt(dx.^2 + dy.^2); 
    
    % Caclulate costs
    % We multiply kappa squared by ds_race to make it grid-independent
    cost_curvature = sum((kappa.^2) .* ds_race); 
    cost_length    = sum(ds_race);               % Total distance traveled
    
    cost = 10*(cost_curvature + (weight_length * cost_length)); % 1000 scaling to reduce function operations
    
end