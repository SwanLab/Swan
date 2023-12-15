% func returns the solution of system with finite difference matrix,
% system is solved by Tridiagonal Matrix Algorithm for pseudo tridiagonal matrix
function v = TriDiagMatrixAlgorithm(obj, g)
    N = size(g, 2);
    % A_i v_{i-1} + B_i v_i + C_i v_{i+1} = g_i
    A = zeros(N, N, N);
    B = zeros(N, N, N);
    B(:, :, 1) = eye(N);
    B(:, :, N) = eye(N);
    C = zeros(N, N, N);
    for i = 2:N-1
        A(:, :, i) = diag([0, obj.d(:, i-1)', 0]);
        B(:, :, i) = diag([1, obj.a(:, i-1)', 1]) +...
            diag([0, obj.e(:, i-1)'], 1) + diag([obj.c(:, i-1)', 0], -1);
        C(:, :, i) = diag([0, obj.b(:, i-1)', 0]);    
    end
    
    alpha = zeros(N, N, N);
    alpha(:, :, 1) = -inv(B(:, :, 1))*C(:, :, 1);
    beta = zeros(N);
    beta(:, 1) = B(:, :, 1)\g(:, 1);
    for i = 2:(N-1)
        D_i = A(:, :, i)*alpha(:, :, i-1) + B(:, :, i);
        alpha(:, :, i) = -inv(D_i)*C(:, :, i);
        beta(:, i) = D_i\(g(:, i) - A(:, :, i)*beta(:, i-1));
    end
    alpha(:, :, N) = zeros(N);
    beta(:, N) = (A(:, :, N)*alpha(:, :, N-1) + B(:, :, N))\(g(:, N) - A(:, :, N)*beta(:, N-1));
    
    v = zeros(N);
    v(:, N) = beta(:, N);
    for i = (N-1):-1:1
        v(:, i) = alpha(:, :, i)*v(:, i+1) + beta(:, i);
    end    
end