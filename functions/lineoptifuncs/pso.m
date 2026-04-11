function [global_solution, gs_eval, gs_history, gs_eval_history] = pso(cost_func, n_dimensions, lb, ub, n_particles, n_iterations, w, cp, cg, verbose)
    % Particle Swarm Optimization - MATLAB Vectorized Version
    
    % Ensure bounds are row vectors for matrix operations
    lb = lb(:)';
    ub = ub(:)';
    v_max = (ub - lb); % Maximum velocity magnitude
    
    % Initialize Swarm Matrices
    % Each row is a particle, each column is a dimension
    positions = lb + (ub - lb) .* rand(n_particles, n_dimensions);
    velocities = -v_max + 2 .* v_max .* rand(n_particles, n_dimensions);
    
    pbest_positions = positions;
    pbest_evals = inf(n_particles, 1);
    
    global_solution = zeros(1, n_dimensions);
    gs_eval = inf;
    
    gs_history = zeros(n_iterations, n_dimensions);
    gs_eval_history = zeros(n_iterations, 1);
    
    if verbose
        fprintf('\n------------------ PARAMETERS -----------------\n');
        fprintf('Number of dimensions: %d\n', n_dimensions);
        fprintf('Number of iterations: %d\n', n_iterations);
        fprintf('Number of particles: %d\n', n_particles);
        fprintf('w: %.2f\tcp: %.2f\tcg: %.2f\n\n', w, cp, cg);
        fprintf('----------------- OPTIMIZATION ----------------\n');
        fprintf('Population initialization...\n');
    end
    
    % Initial Evaluation
    for i = 1:n_particles
        % cost_func expects a column vector, so we transpose
        eval = cost_func(positions(i, :)'); 
        pbest_evals(i) = eval;
        
        if eval < gs_eval
            gs_eval = eval;
            global_solution = positions(i, :);
        end
    end
    
    if verbose
        fprintf('Start of optimization...\n');
    end
    
    tic; % Start timer
    
    % Main PSO Loop
    for k = 1:n_iterations
        % Generate random components for the whole swarm at once
        rp = rand(n_particles, n_dimensions);
        rg = rand(n_particles, n_dimensions);
        
        % Update Velocities (Vectorized)
        velocities = w .* velocities + ...
                     cp .* rp .* (pbest_positions - positions) + ...
                     cg .* rg .* (repmat(global_solution, n_particles, 1) - positions);
                 
        % Clamp Velocities
        velocities = max(min(velocities, v_max), -v_max);
        
        % Update Positions
        positions = positions + velocities;
        
        % Clamp Positions to Lower and Upper Bounds
        positions = max(min(positions, ub), lb);
        
        % Evaluate Swarm
        for i = 1:n_particles
            % Evaluate cost (transpose to column vector)
            p_eval = cost_func(positions(i, :)');
            
            % Update Personal Best
            if p_eval < pbest_evals(i)
                pbest_evals(i) = p_eval;
                pbest_positions(i, :) = positions(i, :);
                
                % Update Global Best
                if p_eval < gs_eval
                    gs_eval = p_eval;
                    global_solution = positions(i, :);
                end
            end
        end
        
        % Store History
        gs_eval_history(k) = gs_eval;
        gs_history(k, :) = global_solution;
        
        if verbose && mod(k, max(1, floor(n_iterations/10))) == 0
            fprintf('Iteration %d/%d - Best Cost: %.5f\n', k, n_iterations, gs_eval);
        end
    end
    
    elapsed_time = toc; % End timer
    
    % Transpose global_solution back to column vector for consistency with your code
    global_solution = global_solution'; 
    
    if verbose
        fprintf('\nEnd of optimization...\n\n');
        fprintf('------------------- RESULTS -------------------\n');
        fprintf('Optimization elapsed time: %.2f s\n', elapsed_time);
        fprintf('Solution evaluation: %.5f\n', gs_eval);
    end
end