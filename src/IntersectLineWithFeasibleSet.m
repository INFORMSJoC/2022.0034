function [ solution ] = IntersectLineWithFeasibleSet(feasible_set, x_bar, x_tilde)

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
        
        %% Find where line through x_tilde and global min pierces polyhedron
        lowerbound = -inf;
        upperbound = inf;
        
        for i = 1 : size(D, 1)
            if D(i, :) * (x_tilde - x_bar) > 1e-6
                upperbound = min(upperbound, (d(i) - D(i,:)*x_bar) / (D(i,:) * (x_tilde - x_bar)));
            elseif D(i, :) * (x_tilde - x_bar) < -1e-6
                lowerbound = max(lowerbound, (d(i) - D(i,:)*x_bar) / (D(i,:) * (x_tilde - x_bar)));
            else
                if d(i) - D(i,:)*x_bar < -1e-6
                    break; % Line does not pierce polyhedron
                end
            end
        end
        
        % Find farthest of the two intersection points
        if lowerbound - upperbound > 1e-6 && (abs(lowerbound) > 1e-3 || abs(upperbound) > 1e-3)
            solution = zeros(n, 1); % Line does not pierce polyhedron
        elseif abs(upperbound) > abs(lowerbound)
            solution = upperbound * x_tilde + ( 1 - upperbound ) * x_bar;
        else
            solution = lowerbound * x_tilde + ( 1 - lowerbound ) * x_bar;
        end
        
    case '1-norm'
        
        solution = x_tilde * feasible_set.rho / norm(x_tilde, 1);
        
    case 'Intersection of Ellipsoids'
        lowerbound = -inf;
        upperbound = inf;
        
        for i = 1 : feasible_set.m
            a_coef = x_tilde'*feasible_set.A{i}*x_tilde - 2 * x_tilde'*feasible_set.A{i}*x_bar + x_bar'*feasible_set.A{i}*x_bar;
            b_coef = 2 * x_tilde'*feasible_set.A{i}*x_bar - 2 * x_bar'*feasible_set.A{i}*x_bar + feasible_set.b{i}'*x_tilde - feasible_set.b{i}'*x_bar;
            c_coef = x_bar'*feasible_set.A{i}*x_bar + feasible_set.b{i}'*x_bar + feasible_set.c{i};
            
            deter = b_coef^2 - 4 * a_coef * c_coef;
            
            if deter > 0
                lowerbound = max(lowerbound, (- b_coef - sqrt(deter)) / (2 * a_coef));
                upperbound = min(upperbound, (- b_coef + sqrt(deter)) / (2 * a_coef));
            elseif deter == 0
                lowerbound = - b_coef / (2 * a_coef);
                upperbound = - b_coef / (2 * a_coef);
            else
                break; % Line does not pierce
            end
        end
        
        % Find farthest of the two intersection points
        if upperbound < lowerbound
            solution = zeros(size(x_tilde)); % Line does not pierce
        elseif abs(upperbound) > abs(lowerbound)
            solution = upperbound * x_tilde + ( 1 - upperbound ) * x_bar;
        else
            solution = lowerbound * x_tilde + ( 1 - lowerbound ) * x_bar;
        end
        
    case 'box'
        n = feasible_set.n;
        D = [-eye(n); eye(n)];
        d = [-feasible_set.lb; feasible_set.ub];
        
        %% Find where line through x_tilde and global min pierces polyhedron
        lowerbound = -inf;
        upperbound = inf;
        
        for i = 1 : size(D, 1)
            if D(i, :) * (x_tilde - x_bar) > 0
                upperbound = min(upperbound, (d(i) - D(i,:)*x_bar) / (D(i,:) * (x_tilde - x_bar)));
            elseif D(i, :) * (x_tilde - x_bar) < 0
                lowerbound = max(lowerbound, (d(i) - D(i,:)*x_bar) / (D(i,:) * (x_tilde - x_bar)));
            else
                if d(i) - D(i,:)*x_bar < 0
                    break; % Line does not pierce
                end
            end
        end
        
        % Find farthest of the two intersection points
        if upperbound < lowerbound
            solution = zeros(size(x_tilde)); % Line does not pierce
        elseif abs(upperbound) > abs(lowerbound)
            solution = upperbound * x_tilde + ( 1 - upperbound ) * x_bar;
        else
            solution = lowerbound * x_tilde + ( 1 - lowerbound ) * x_bar;
        end
        
    case 'Aras 4.1'
        if abs( norm(x_tilde - feasible_set.a, 2) - feasible_set.rho) < 1e-5
            solution = x_tilde;
        else
            keyboard
        end
        
    otherwise
        error('Line method is not defined for this type of feasible set.')
        
end

end
