function Z = apply_PA(PA, U, Z)
% Apply piecewise-linear Jacobian to columns Z (same size as U)
% U: original unconstrained columns (m x k)
% Z: derivative columns to be mapped (m x k), overwritten in place

    [m,k] = size(U);
    % Determine active sets: A_c = {i : U_i > theta_c}
    A = U > PA.theta;                        % (m x k) logical
    % Sum of Z over active sets (per column)
    sZ = sum(Z .* A, 1);                     % (1 x k)
    inv_rho = 1 ./ max(PA.rho, 1);           % guard rho=0 (won’t happen if V>0)
    % Broadcast subtract average over active set
    Z = Z .* A - (ones(m,1) * (inv_rho .* sZ));  % zero outside A automatically wrong?
    % Zero outside A explicitly to be safe:
    Z(~A) = 0;
end
