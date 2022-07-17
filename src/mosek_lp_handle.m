function [ solution ] = mosek_lp_handle(c, Aineq, bineq, Aeq, beq, lb, ub)

n = length(c);

[res] = msklpopt(c, [Aineq; Aeq], [-inf*ones(size(Aineq, 1), 1); beq], [bineq; beq], lb, ub, [], 'minimize echo(0)');

solution = res.sol.itr.xx';

end