function [ final_x, full_solutions, computation_times, groups ] = CoMax( objective, feasible_set, comax_settings )
% CoMax - Algorithm to maximize a convex function over a bounded convex set.
%
% INPUT ARGUMENTS
% -----------------------------------------------------------------------------
% objective     - Structure that describes the objective through these fields:
% f         - Function handle that evaluates the objective function at a given
% input x.
% grad      - Function handle that evaluates the gradient of the objective
% function at a given input x.
% -----------------------------------------------------------------------------
% feasible_set  - Structure that describes the feasible set through these
% fields:
% max_linear        - Function that maximizes a linear function with the first
% argument as coefficients over the feasible set and returns the optimal solution.
% random_boundary   - Function that returns a random point on the boundary
% of the feasible set.
% n                 - Dimension of the feasible set.
% -----------------------------------------------------------------------------
% comax_settings          - Structure that allows for the following fields:
% max_iterations    - The maximum number of iterations.
% stop_threshold    - The maximum 2-norm difference between two subsequent
% solutions that terminates the algorithm.
% N                 - Number of random starting points from which the algorithm
% should be run.
% -----------------------------------------------------------------------------

%% Process Settings
if nargin < 2
    error("Not enough input arguments supplied.")
elseif nargin < 3
    comax_settings = struct();
end

if ~isfield(comax_settings, 'max_iterations')
    comax_settings.max_iterations = 50;
end
if ~isfield(comax_settings, 'stop_threshold')
    comax_settings.stop_threshold = 1e-5;
end
if ~isfield(comax_settings, 'N')
    comax_settings.N = 100;
end

%% Process Input
% Define step function
if isfield(feasible_set, 'max_linear') && isfield(objective, 'grad')
    step = @ (x) feasible_set.max_linear(objective.grad(x));
else
    error("You need to provide a function that returns the maximum of a linear function over the feasible set as well as the gradient of the objective.")
end

n = feasible_set.n;

if n <= 25
    startMethods = {'Box', 'Analytic Center', 'Random', 'MVIE'};
    N = 5 + 10 + comax_settings.N + 6;
    groups = [0, 5, 15, 15 + comax_settings.N, 15 + comax_settings.N + 6];
elseif n <= 100
    startMethods = {'Box', 'Analytic Center', 'Random'};
    N = 5 + 10 + comax_settings.N;
    groups = [0, 5, 15, 15 + comax_settings.N];
else
    startMethods = {'Box', 'Analytic Center'};
    N = 5 + 10;
    groups = [0, 5, 15];
end

%% Initialize Containers

x = zeros(n, comax_settings.max_iterations, N);

final_x = zeros(n, N);

final_it = zeros(N, 1);

computation_times = zeros(N, 1);

%% Run algorithm for all starting methods
for methodIndex = 1 : length(startMethods)
    tic
    
    %% Phase 1
    switch startMethods{methodIndex}
        case 'Box'
            [x(:, 1, 1), x(:, 1, 2), x(:, 1, 3), x(:, 1, 4), x(:, 1, 5)] = ...
                BoundingBox(objective, feasible_set);
            
            startIndex = 1; endIndex = 5;
            
        case 'Analytic Center'
            [x(:, 1, 6), x(:, 1, 7), x(:, 1, 8), x(:, 1, 9), x(:, 1, 10), x(:, 1, 11), ...
                x(:, 1, 12), x(:, 1, 13), x(:, 1, 14), x(:, 1, 15)] = ...
                AnalyticCenterStart(objective, feasible_set);
            
            startIndex = 6; endIndex = 15;
            
        case 'Random'
            for i = 1 : comax_settings.N
                x(:, 1, 15 + i) = feasible_set.random_boundary();
            end
            
            startIndex = 16; endIndex = 15 + comax_settings.N;
            
        case 'MVIE'
            [x(:, 1, 15 + comax_settings.N + 1), x(:, 1, 15 + comax_settings.N + 2), x(:, 1, 15 + comax_settings.N + 3), ...
                x(:, 1, 15 + comax_settings.N + 4), x(:, 1, 15 + comax_settings.N + 5), x(:, 1, 15 + comax_settings.N + 6)] = ...
                MVIE(objective, feasible_set);
            
            startIndex = 15 + comax_settings.N + 1; endIndex = 15 + comax_settings.N + 6;
            
        otherwise
            error('This type of initialization is unknown.')
    end
    
    %% Phase 2
    for i = startIndex : endIndex
        % First iteration
        k = 2;
        
        x(:, k, i) = step(x(:, k-1, i));
        
        % Main loop
        k = 3;
        while norm(x(:, k-1, i) - x(:, k-2, i), 2) > comax_settings.stop_threshold && k <= comax_settings.max_iterations
            x(:, k, i) = step(x(:, k-1, i));
            
            k = k + 1;
        end
        
        % Determine end iteration
        k = k - 1;
        while any(isnan(x(:, k, i)))
            k = k - 1;
        end
        
        final_it(i) = k;
    
        final_x(:, i) = x(:, final_it(i), i);
    end
    
    computation_times(methodIndex) = toc;
    
end

full_solutions = cell(N, 1);

for i = 1 : N
    full_solutions{i} = x(:, 1:final_it(i), i);
end

end