function [ line_inscribed, direction_inscribed, gradient_inscribed, ...
    line_inscribed_cons, direction_inscribed_cons, gradient_inscribed_cons, ...
    direction_circumscribing, gradient_circumscribing, ...
    direction_circumscribing_cons, gradient_circumscribing_cons] = ...
    AnalyticCenterStart( objective, feasible_set )

analytic_center_settings = sdpsettings('solver', 'mosek', 'verbose', 0);

%% Find Analytic Center

switch feasible_set.type
    case 'polyhedron'
        n = feasible_set.n;
        
        D = feasible_set.D;
        d = feasible_set.d;
        
        I = eye(n);
        
        if ~ isempty(feasible_set.lb)
            D = [D; - I(~isinf(feasible_set.lb), :)];
            d = [d; - feasible_set.lb(~isinf(feasible_set.lb))];
        end
        if ~ isempty(feasible_set.ub)
            D = [D; I(~isinf(feasible_set.ub), :)];
            d = [d; feasible_set.ub(~isinf(feasible_set.ub))];
        end
        
        m = size(D, 1);
        
        x = sdpvar(n, 1);
        
        result = optimize([], -sum(log(d - D * x)), analytic_center_settings);
        
        center = double(x);
        
    case '1-norm'
        
        center = zeros(feasible_set.n, 1);
        
    case 'Intersection of Ellipsoids'
        n = feasible_set.n;
        m = feasible_set.m;
        
        %% Norm based formulation for AC
        
        % Order of variables x, y, t

        [~, res_symb] = mosekopt('symbcon');

        % Objective
        prob.c      = [zeros(1, n), zeros(1, m), ones(1, m)];
        
        % Linear Inequalities
        prob.a      = sparse(zeros(0, n + m + m));    
        
        % Exp Cone Inequalities
        prob.cones  = repmat([res_symb.symbcon.MSK_CT_PEXP 3], [1 m]);
        prob.f      = sparse([1:3:3*m, 3:3:3*m], [n + (1 : m), n + m + (1 : m)], [ones(1, m), ones(1, m)]);
        prob.g      = repmat([0; 1; 0], [m 1]);
        
        % Original Feasible Set Inequalities
        prob.cones  = [prob.cones, repmat([res_symb.symbcon.MSK_CT_QUAD n+2], [1 m])];
        for i = 1 : m
            prob.f  = [prob.f; sparse([
                zeros(1, n), zeros(1, i-1), -1/2, zeros(1, m-i), zeros(1, m);
                feasible_set.Q{i}, zeros(n, m), zeros(n, m);
                zeros(1, n), zeros(1, i-1), 1/2, zeros(1, m-i), zeros(1, m)])];
            prob.g  = [prob.g; 1; feasible_set.q{i}; 0];
        end
        
        [r, res] = mosekopt('maximize echo(0)', prob);

        center = res.sol.itr.xx(1:n);
        
    case 'box'
        center = zeros(feasible_set.n, 1);
        
    case 'Aras 4.1'
        center = feasible_set.a;
end


%% Define Hessian
switch feasible_set.type
    case 'polyhedron'
        Hess = D'* diag( (1 ./ (d - D * center) ).^2 ) * D;
        
        inscr_radius = 1;
%         circum_radius = m * (m-1);
        circum_radius = (1 + 2/sqrt(m)) * m;
        
    case '1-norm'
        Hess = eye(feasible_set.n);
        
        inscr_radius = 1/feasible_set.n;
        circum_radius = 1;
        
    case 'Intersection of Ellipsoids'
        Hess = zeros(n);
        
        for i = 1 : m
            Hess = Hess - (2 * feasible_set.A{i} * (center' * feasible_set.A{i} * center + feasible_set.b{i}'*center + feasible_set.c{i}) - ...
                ( (2 * feasible_set.A{i} * center + feasible_set.b{i} ) * (2 * feasible_set.A{i} * center + feasible_set.b{i} )' )) / ...
                ( (center' * feasible_set.A{i} * center + feasible_set.b{i}'*center + feasible_set.c{i}) ^2 );
        end
        
        inscr_radius = 1;
        circum_radius = (1 + 2/sqrt(m)) * m;
        
    case 'box'
        Hess = eye(feasible_set.n);
        inscr_radius = 1;
        circum_radius = sqrt(feasible_set.n);
        
    case 'Aras 4.1'
        Hess = eye(feasible_set.n);
        inscr_radius = feasible_set.rho;
        circum_radius = feasible_set.rho;
        
end

%% Inscribed Ellipsoid (Unconstrained)
x_tilde = FarthestOnEllipsoid(objective.global_min, objective, 2*Hess, (-2*center'*Hess)', center'*Hess*center - inscr_radius);

line_inscribed = IntersectLineWithFeasibleSet(feasible_set, objective.global_min, x_tilde);
direction_inscribed = feasible_set.max_linear(x_tilde - objective.global_min);
gradient_inscribed = feasible_set.max_linear(objective.grad(x_tilde));


%% Inscribed Ellipsoid (Constrained)
x_tilde = FarthestOnEllipsoid(objective.x_bar_cons, objective, 2*Hess, (-2*center'*Hess)', center'*Hess*center - inscr_radius);

line_inscribed_cons = IntersectLineWithFeasibleSet(feasible_set, objective.x_bar_cons, x_tilde);
direction_inscribed_cons = feasible_set.max_linear(x_tilde - objective.x_bar_cons);
gradient_inscribed_cons = feasible_set.max_linear(objective.grad(x_tilde));


%% Circumscribing Ellipsoid (Unconstrained)
x_tilde = FarthestOnEllipsoid(objective.global_min, objective, 2*Hess, (-2*center'*Hess)', center'*Hess*center - circum_radius);

direction_circumscribing = feasible_set.max_linear(x_tilde - objective.global_min);
gradient_circumscribing = feasible_set.max_linear(objective.grad(x_tilde));

%% Circumscribing Ellipsoid (Constrained)

x_tilde = FarthestOnEllipsoid(objective.x_bar_cons, objective, 2*Hess, (-2*center'*Hess)', center'*Hess*center - circum_radius);

direction_circumscribing_cons = feasible_set.max_linear(x_tilde - objective.x_bar_cons);
gradient_circumscribing_cons = feasible_set.max_linear(objective.grad(x_tilde));

end