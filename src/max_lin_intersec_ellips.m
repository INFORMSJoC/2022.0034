function [ solution ] = max_lin_intersec_ellips(y, feasible_set)

feasible_set.mosek_prob.c = y';

[r, res] = mosekopt('maximize echo(0)', feasible_set.mosek_prob);

solution = res.sol.itr.xx;

feasible_set.mosek_prob.c = [];


end

