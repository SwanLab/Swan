function [w_after, feasible] = nn_update_lsqlin(Uloc, b, w_before)
    % Minimize ||w - w_before||^2 s.t. Uloc' * w = b, w >= 0
    n = length(w_before);
    C = speye(n);          % objective: ||C*w - d||^2 with d = w_before
    d = w_before;
    Aeq = Uloc'; beq = b;
    lb = zeros(n,1);       % nonnegativity
    ub = [];               % no upper bounds

    opts = optimoptions('lsqlin','Display','off',...
                        'Algorithm','interior-point'); % robust with eq+lb
    [w_after, ~, residual_norm, exitflag] = lsqlin(C, d, [], [], Aeq, beq, lb, ub, [], opts);
%     
%     residual = norm(Uloc' * w_after - b);
%     tol = 1e-10 * norm(b);
%     feasible = (residual <= tol) && all(w_after >= -1e-12) && (exitflag > 0);
if exitflag ==1 
    feasible = true; 
else
    feasible = false; 
end
end
