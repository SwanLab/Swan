function [W, PA] = project_simplex_cols(U, V)
%PROJECT_SIMPLEX_COLS  Euclidean projection, column-wise, onto {w>=0, sum w = V}
%   U: (m x k) unconstrained columns
%   V: scalar target sum
%   W: (m x k) projected
%   PA: struct enabling fast Jacobian action on derivatives (see below)

    [m,k] = size(U);
    % Sort each column descending
    [Us, idx] = sort(U, 1, 'descend');          % (m x k)
    css = cumsum(Us,1);                          % (m x k)
    j = (1:m)';                                  % (m x 1)
    T = (css - V) ./ j;                          % (m x k)
    cond = Us - T > 0;                           % (m x k) logical
    % rho: last true per column
    rho = sum(cond,1);                           % (1 x k)
    % theta: pick T(rho,c) for each column
    lin = sub2ind([m,k], rho, 1:k);
    theta = T(lin);                               % (1 x k)
    % Project
    W = max(U - theta, 0);

    % --- Compact Jacobian info for derivative mapping (optional) ---
    % Active set sizes
    s = rho;                                     % (1 x k), |A| per column
    % For each column c, the Jacobian P_A acts as:
    %   (P_A * v)_i = v_i - (1/s_c) * sum_{j in A_c} v_j   if i in A_c
    %                  0                                    otherwise
    % We encode what’s needed to apply this without forming P_A:
    PA.theta = theta;     % threshold
    PA.rho   = rho;       % |A|
    PA.idx   = idx;       % sorting permutation
end
