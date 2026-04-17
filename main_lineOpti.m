%% MAIN

clc;
clear;
close all;
echo off; % Used to avoid unwanted warnings and other console stuff

addpath("functions\")
addpath("functions\lineoptifuncs\")

%% MAIN
clc;
clear;
close all;
echo off; % Used to avoid unwanted warnings and other console stuff

addpath("functions\")

% Run Params
trackplot=false;

par = carParams();
n_var = 200;          % Number of Design Variables for Interpolation
car_margin = 0.5;    % Car half-width margin (e.g., 1 meter wide car = 0.5m margin)

% Generate Track
track = genTrack(trackplot);

% Set up for optimizer
lineopti.s_full = [0; cumsum(track.vecmag(1:end-1))];        % Cumulative distance
lineopti.s_ctrl = linspace(0, lineopti.s_full(end), n_var)';

lineopti.w_left_ctrl = interp1(lineopti.s_full, track.w(:,1), lineopti.s_ctrl);
lineopti.w_right_ctrl = interp1(lineopti.s_full, track.w(:,2), lineopti.s_ctrl);

% fmincon lower and upper bounds (alpha is positive to the left)

lb = -lineopti.w_right_ctrl + car_margin; 
ub =  lineopti.w_left_ctrl  - car_margin;

% Initial guess (start exactly on the centerline, so alpha = 0)
lineopti.alpha_guess = zeros(n_var, 1);

% figure;
% plot(lineopti.s_full,track.w(:,1))
% hold on
% plot(lineopti.s_ctrl,lineopti.w_left_ctrl)

%%% OPTIMIZATION ROUTINE

% Define fmincon options
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'iter', ...
    'MaxFunctionEvaluations', 2000000, ... 
    'MaxIterations', 500, ...           
    'StepTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6);

weight_length = 0.001;
objectiveFcn = @(alpha) calcCurvatureCost(alpha, lineopti.s_ctrl, lineopti.s_full, track, weight_length);

% Equality constraints
% First and last point need to have the same position and direction
Aeq = zeros(2, n_var);
beq = zeros(2, 1); % The right side of the equations (both equal 0)

Aeq(1, 1)   = 1;
Aeq(1, end) = -1;

Aeq(2, 1)     = -1;
Aeq(2, 2)     =  1;
Aeq(2, end-1) =  1;
Aeq(2, end)   = -1;

% x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
lineopti.alpha_opt = fmincon(objectiveFcn, lineopti.alpha_guess, [], [], Aeq, beq, lb, ub, [], options);

%%% Generate optimized line and plot result
lineopti.alpha_opti_full = makima(lineopti.s_ctrl, lineopti.alpha_opt, lineopti.s_full);

bdeg = 3;
bknots = augknt(lineopti.s_ctrl,bdeg+1);
b_spline_curve = spmak(bknots,lineopti.alpha_opt');
lineopti.alpha_opti_full = fnval(b_spline_curve, lineopti.s_full);

lineopti.optimized = track.m+track.vecleft./track.vecmag.*lineopti.alpha_opti_full;

figure;
plot(track.m(:,1),track.m(:,2),lineopti.optimized(:,1),lineopti.optimized(:,2),track.l(:,1),track.l(:,2),track.r(:,1),track.r(:,2))
legend("Centerline","Optimized Line","Left Bound","Right Bound")

%%% Calculate Line Parameters

lineopti.dx = gradient(lineopti.optimized(:,1));
lineopti.dy = gradient(lineopti.optimized(:,2));
lineopti.ddx = gradient(lineopti.dx);
lineopti.ddy = gradient(lineopti.dy);
lineopti.kappasquare = sum(((lineopti.dx .* lineopti.ddy - lineopti.dy .* lineopti.ddx) ./ ((lineopti.dx.^2 + lineopti.dy.^2).^(3/2))).^2);
lineopti.length = sum(sqrt(lineopti.dx.^2 + lineopti.dy.^2));

fprintf("Track Midline Square Curvature: %e\n", track.kappasquare)
fprintf("Track Midline Length [m]: %e\n", track.length)
fprintf("Optimized Line Square Curvature: %e\n", lineopti.kappasquare)
fprintf("Line Length [m]: %e\n", lineopti.length)

[laptime,vprof] = calcTimeAndVelocity(lineopti.alpha_opt, lineopti.s_ctrl, lineopti.s_full, track, par, 'bspline');
fprintf("Lap Time [s]: %e\n", laptime)

figure;
plot(lineopti.s_full,vprof)
% test2 = calcLapTimeCost(lineopti.alpha_opt, lineopti.s_ctrl, lineopti.s_full, track, par)
% 
% figure;
% plot(test1)
% hold on
% plot(test2)
% 
% a = mean(test1)
% b = mean(test2)