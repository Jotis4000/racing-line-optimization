%% MAIN
clc;
clear;
close all;
echo off; 
addpath("functions\")
addpath("functions\lineoptifuncs\")

% Run Params
trackplot = false;
par = carParams();
splineType = 'makima'; % Keeping 'makima' based on your previous bound-overshoot fix

% -------------------------------------------------------------------------
% MULTI-GRID SETUP
n_var_pso  = 120;   % Coarse Grid: Fast global search
n_var_fmin = 350;  % Fine Grid: High-resolution local smoothing
% -------------------------------------------------------------------------
car_margin = 0.5;  

% Generate Track
track = genTrack(trackplot);

% Set up baseline distance
lineopti.s_full = [0; cumsum(track.vecmag(1:end-1))]; 

%%% =======================================================
%%% ADAPTIVE NODE SPACING (Creating TWO Grids)
%%% =======================================================
corner_weight = 20;    
smoothing_window = 30; 

% 1. Coarse Grid (For PSO)
lineopti.s_ctrl_pso = genAdaptiveNodes(lineopti.s_full, track.m, n_var_pso, corner_weight, smoothing_window);
lineopti.w_left_pso  = interp1(lineopti.s_full, track.w(:,1), lineopti.s_ctrl_pso, 'spline');
lineopti.w_right_pso = interp1(lineopti.s_full, track.w(:,2), lineopti.s_ctrl_pso, 'spline');

lb_pso = -lineopti.w_right_pso + car_margin; 
ub_pso =  lineopti.w_left_pso  - car_margin;

Aeq_pso = zeros(2, n_var_pso);
Aeq_pso(1, 1) = 1; Aeq_pso(1, end) = -1;
Aeq_pso(2, 1) = -1; Aeq_pso(2, 2) = 1; Aeq_pso(2, end-1) = 1; Aeq_pso(2, end) = -1;

% 2. Fine Grid (For fmincon)
lineopti.s_ctrl_fmin = genAdaptiveNodes(lineopti.s_full, track.m, n_var_fmin, corner_weight, smoothing_window);
lineopti.w_left_fmin  = interp1(lineopti.s_full, track.w(:,1), lineopti.s_ctrl_fmin, 'spline');
lineopti.w_right_fmin = interp1(lineopti.s_full, track.w(:,2), lineopti.s_ctrl_fmin, 'spline');

lb_fmin = -lineopti.w_right_fmin + car_margin; 
ub_fmin =  lineopti.w_left_fmin  - car_margin;

Aeq_fmin = zeros(2, n_var_fmin);
Aeq_fmin(1, 1) = 1; Aeq_fmin(1, end) = -1;
Aeq_fmin(2, 1) = -1; Aeq_fmin(2, 2) = 1; Aeq_fmin(2, end-1) = 1; Aeq_fmin(2, end) = -1;

beq = zeros(2, 1); 
%%% =======================================================


%%% =======================================================
%%% THE HYBRID SOLVER: Cascaded Optimization
%%% =======================================================

%%% PHASE 1: GLOBAL SEARCH (COARSE GRID)
fprintf('\n--- Starting Phase 1: PSO (Global Search on %d nodes) ---\n', n_var_pso);
n_particles = 200;     
n_iterations = 1000;   
w = 0.7; cp = 1.5; cg = 1.5; 
verbose = true;       

penalty_weight = 1e5;
% Using your new friction cost function for the optimizer
objectiveFcnPSO = @(alpha) calcLapTimeCostDetailFriction(alpha, lineopti.s_ctrl_pso, lineopti.s_full, track, par, splineType) ...
                           + penalty_weight * sum((Aeq_pso * alpha - beq).^2);

[alpha_pso_rough, best_cost_pso, gs_hist, eval_hist] = pso(...
    objectiveFcnPSO, n_var_pso, lb_pso, ub_pso, n_particles, n_iterations, w, cp, cg, verbose);

fprintf('\nPhase 1 Complete. Bridging grids...\n');

%%% THE BRIDGE: Upscale 50 nodes to 150 nodes
alpha_fmin_guess = makima(lineopti.s_ctrl_pso, alpha_pso_rough, lineopti.s_ctrl_fmin);
alpha_fmin_guess = max(min(alpha_fmin_guess, ub_fmin), lb_fmin); % Safety clamp

%%% PHASE 2: LOCAL POLISH (FINE GRID)
fprintf('\n--- Starting Phase 2: fmincon (Local Smoothing on %d nodes) ---\n', n_var_fmin);

objectiveFcnExact = @(alpha) calcLapTimeCostDetailFriction(alpha, lineopti.s_ctrl_fmin, lineopti.s_full, track, par, splineType);

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...       
    'Display', 'iter', ...
    'MaxFunctionEvaluations', 20000000, ...
    'MaxIterations', 2000, ...
    'StepTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6);

lineopti.alpha_opt = fmincon(objectiveFcnExact, alpha_fmin_guess, [], [], Aeq_fmin, beq, lb_fmin, ub_fmin, [], options);

fprintf('\n--- Hybrid Optimization Complete ---\n');


%%% =======================================================
%%% POST-PROCESSING & PLOTTING (Updated to your format)
%%% =======================================================

%%% Generate optimized line
if isequal(splineType,'makima')
    lineopti.alpha_opti_full = makima(lineopti.s_ctrl_fmin, lineopti.alpha_opt, lineopti.s_full);
elseif isequal(splineType,'bspline')
    bdeg = 3;
    bknots = augknt(lineopti.s_ctrl_fmin, bdeg+1);
    b_spline_curve = spmak(bknots, lineopti.alpha_opt');
    lineopti.alpha_opti_full = fnval(b_spline_curve, lineopti.s_full);
end

lineopti.optimized = track.m + track.vecleft ./ track.vecmag .* lineopti.alpha_opti_full;

% Your requested Plot format
figure;
plot(track.m(:,1), track.m(:,2), lineopti.optimized(:,1), lineopti.optimized(:,2), track.l(:,1), track.l(:,2), track.r(:,1), track.r(:,2))
legend("Centerline", "Optimized Line", "Left Bound", "Right Bound")

%%% Calculate Line Parameters
lineopti.dx = gradient(lineopti.optimized(:,1));
lineopti.dy = gradient(lineopti.optimized(:,2));
lineopti.ddx = gradient(lineopti.dx);
lineopti.ddy = gradient(lineopti.dy);
lineopti.kappasquare = sum(((lineopti.dx .* lineopti.ddy - lineopti.dy .* lineopti.ddx) ./ ((lineopti.dx.^2 + lineopti.dy.^2).^(3/2))).^2);
lineopti.ds = sqrt(lineopti.dx.^2 + lineopti.dy.^2);
lineopti.length = sum(lineopti.ds);

fprintf("Track Midline Square Curvature: %e\n", track.kappasquare)
fprintf("Track Midline Length [m]: %e\n", track.length)
fprintf("Optimized Line Square Curvature: %e\n", lineopti.kappasquare)
fprintf("Line Length [m]: %e\n", lineopti.length)

% Calculate final actual laptime using your dedicated calcTimeAndVelocity function
% Notice we pass s_ctrl_fmin here because lineopti.alpha_opt is mapped to the fine grid!
[laptime, vprof] = calcTimeAndVelocity(lineopti.alpha_opt, lineopti.s_ctrl_fmin, lineopti.s_full, track, par, splineType);

fprintf("Lap Time [s]: %e\n", laptime)

% Your requested Velocity Plot format
figure;
plot(lineopti.s_full, vprof)