%% MAIN
clc;
clear;
close all;
echo off; % Used to avoid unwanted warnings and other console stuff

addpath("functions\")
addpath("functions\lineoptifuncs\")

% Run Params
trackplot=false;

par = carParams();
n_var = 200;          % Number of Design Variables for Interpolation
car_margin = 0.5;    % Car half-width margin (e.g., 1 meter wide car = 0.5m margin)

guess = 'geom'; % 'center', 'geom'
splineType = 'makima'; % 'makima', 'bspline'

% Generate Track
track = genTrack(trackplot);

% Set up for optimizer
lineopti.s_full = [0; cumsum(track.vecmag(1:end-1))];        % Cumulative distance
lineopti.s_ctrl = linspace(0, lineopti.s_full(end), n_var)';

corner_weight = 20;    % How strongly corners pull nodes (Multiplier)
smoothing_window = 50; % How wide the apex cluster is
lineopti.s_ctrl = genAdaptiveNodes(lineopti.s_full, track.m, n_var, corner_weight, smoothing_window);

lineopti.w_left_ctrl = interp1(lineopti.s_full, track.w(:,1), lineopti.s_ctrl);
lineopti.w_right_ctrl = interp1(lineopti.s_full, track.w(:,2), lineopti.s_ctrl);

% fmincon lower and upper bounds (alpha is positive to the left)

if isequal(splineType,'makima')
    bmargin = 0;
elseif isequal(splineType,'bspline')
    bmargin = 10;
end
lb = -lineopti.w_right_ctrl + car_margin-bmargin; 
ub =  lineopti.w_left_ctrl  - car_margin+bmargin;

%%% OPTIMIZATION ROUTINE



% Initial guess (start exactly on the centerline, so alpha = 0)
lineopti.alpha_guess = zeros(n_var, 1);

if isequal(guess,'center')
    lineopti.alpha_guess = zeros(n_var, 1);
elseif isequal(guess,'geom')
    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'iter', ...
        'MaxFunctionEvaluations', 2000000, ... 
        'MaxIterations', 100, ...           
        'StepTolerance', 1e-6, ...
        'OptimalityTolerance', 1e-6);

    weight_length = 0;
    objectiveFcn = @(alpha) calcCurvatureCost(alpha, lineopti.s_ctrl, lineopti.s_full, track, weight_length);

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
    lineopti.alpha_guess = lineopti.alpha_opt;
end

random = randn(n_var);
pert = random(:,1);

stdv = 0.5;
lineopti.alpha_guess = lineopti.alpha_guess+0.5*pert;
for i=1:n_var
    if lineopti.alpha_guess(i)<lb(i)
        lineopti.alpha_guess(i)=lb(i)+0.01;
    end
    % lineopti.alpha_guess(i)
    if lineopti.alpha_guess(i)>ub(i)
        lineopti.alpha_guess(i)=ub(i)-0.01;
    end
end

plot(lineopti.s_ctrl,lineopti.alpha_guess-0.5*pert,lineopti.s_ctrl,lineopti.alpha_guess,lineopti.s_ctrl,lb,lineopti.s_ctrl,ub)