function [ start_line, start_direction, start_gradient, ...
    start_line_cons, start_direction_cons, start_gradient_cons ] = ...
    MVIE(objective, feasible_set)

ellipsoid_settings = sdpsettings('solver', 'mosek', 'verbose', 0);

%% First find maximimum volume inscribed ellipsoid

switch feasible_set.type
    case 'polyhedron'
        n = feasible_set.n;
        
        D = feasible_set.D;
        d = feasible_set.d;
        if (~ isempty(feasible_set.lb)) && (~ all(isinf(-feasible_set.lb)))
            D = [D; -eye(n)];
            d = [d; -feasible_set.lb];
        end
        
        if (~ isempty(feasible_set.ub)) && (~ all(isinf(feasible_set.ub)))
            D = [D; eye(n)];
            d = [d; feasible_set.ub];
        end
        
        
        Q = sdpvar(n);
        
        q = sdpvar(n, 1, 'full');
        
        F = [];
        for i = 1 : length(d)
            F = [F; norm(Q * D(i, :)', 2) + D(i, :) * q <= d(i)];
        end
        
        obj = - logdet(Q);
        
        result = optimize(F, obj, ellipsoid_settings);
        
        Q = double(Q);
        q = double(q);
        
    case '1-norm'
        n = feasible_set.n;
        
        Q = feasible_set.rho * 1/n * sqrt(n) * eye(n);
        q = zeros(n, 1);
        
    case 'Intersection of Ellipsoids'
        n = feasible_set.n;
        
        Q = sdpvar(n);
        q = sdpvar(n, 1, 'full');
        
        lambda = sdpvar(feasible_set.m, 1, 'full');
        
        obj = - logdet(Q);
        
        F = [];
        for i = 1 : feasible_set.m
            F = [F;
                [-lambda(i) - feasible_set.c{i} + 1/2 * feasible_set.b{i}' * inv(feasible_set.A{i}) * 1/2*feasible_set.b{i}, zeros(1, n), (q + inv(feasible_set.A{i}) * 1/2*feasible_set.b{i})';
                zeros(n, 1), lambda(i) * eye(n), Q;
                q + inv(feasible_set.A{i}) * 1/2*feasible_set.b{i}, Q, inv(feasible_set.A{i})] >= 0];
        end
        
        result = optimize(F, obj, ellipsoid_settings);
        
        Q = double(Q);
        q = double(q);
        
    case 'box'
        n = feasible_set.n;
        
        Q = eye(n);
        q = zeros(n, 1);
end

yalmip('clear')

%% Unconstrained initialization
% Find point farthest on ellipsoid

x_tilde = FarthestOnEllipsoid(objective.global_min, objective, 2*inv(Q)*inv(Q), (-2*q'*inv(Q)*inv(Q))', q'*inv(Q)*inv(Q)*q - 1);

start_line = IntersectLineWithFeasibleSet(feasible_set, objective.global_min, x_tilde);
start_direction = feasible_set.max_linear(x_tilde - objective.global_min);
start_gradient = feasible_set.max_linear(objective.grad(x_tilde));

%% Constrained initialization

x_tilde = FarthestOnEllipsoid(objective.x_bar_cons, objective, 2*inv(Q)*inv(Q), (-2*q'*inv(Q)*inv(Q))', q'*inv(Q)*inv(Q)*q - 1);

start_line_cons = IntersectLineWithFeasibleSet(feasible_set, objective.x_bar_cons, x_tilde);
start_direction_cons = feasible_set.max_linear(x_tilde - objective.x_bar_cons);
start_gradient_cons = feasible_set.max_linear(objective.grad(x_tilde));

end

