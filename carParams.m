%%% SET CAR PARAMETERS

function par = carParams()

    % Simple Params
    par.Vmax = 80; % m/s
    par.m = 1000; % kg
    
    % Second Order Params - Required for more advanced model
    par.CLA = 6;
    par.accelmax = 2; % g
    par.decelmax = 4; % g
    par.latmax = 3; % g

end