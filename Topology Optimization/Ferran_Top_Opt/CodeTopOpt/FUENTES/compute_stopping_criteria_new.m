function [stopping_criteria_opt] = compute_stopping_criteria_new(theta,incre_gamma,constraints,optimality_tol,constr_tol,algorithm,penalty)
global iter
switch algorithm
    case 'level_set'
        active_constr = penalty > 0;
        stopping_criteria_opt = theta >= optimality_tol || any(abs(constraints(active_constr)) > constr_tol(active_constr));% || iter < 10;

    case 'Projected_gradient'
        active_constr = penalty > 0;
        stopping_criteria_opt = incre_gamma >= optimality_tol || any(abs(constraints(active_constr)) > constr_tol(active_constr));% || iter < 10;
        
end

end