function [ start_direction, start_gradient, ...
    start_line_cons, start_direction_cons, start_gradient_cons ] = ...
    BoundingBox(objective, feasible_set)

minima = zeros(feasible_set.n, 1);
maxima = zeros(feasible_set.n, 1);

I = eye(feasible_set.n);

for i = 1 : feasible_set.n
    x = feasible_set.max_linear(I(:, i));
    maxima(i) = x(i);
    
    x = feasible_set.max_linear(-I(:, i));
    minima(i) = x(i);
end

%% Unconstrained initialization
x_tilde = zeros(feasible_set.n, 1);
for i = 1 : feasible_set.n
    if abs(maxima(i) - objective.global_min(i)) > abs(minima(i) - objective.global_min(i))
        x_tilde(i) = maxima(i);
    else
        x_tilde(i) = minima(i);
    end
end

start_direction = feasible_set.max_linear(x_tilde - objective.global_min);
start_gradient = feasible_set.max_linear(objective.grad(x_tilde));

%% Constrained initialization

x_tilde = zeros(feasible_set.n, 1);
for i = 1 : feasible_set.n
    if abs(maxima(i) - objective.x_bar_cons(i)) > abs(minima(i) - objective.x_bar_cons(i))
        x_tilde(i) = maxima(i);
    else
        x_tilde(i) = minima(i);
    end
end

start_line_cons = intersect_line_with_feasible_set(feasible_set, objective.x_bar_cons, x_tilde);
start_direction_cons = feasible_set.max_linear(x_tilde - objective.x_bar_cons);
start_gradient_cons = feasible_set.max_linear(objective.grad(x_tilde));

end