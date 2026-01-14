function [Q,M] =  OrthPoly_ChatGPT2(xi,max_degree,PLOT_FUNCTIONS,ORTHOGONALITY_L2)
% 1. Define the mesh
%n_elements = 4;
%x_nodes = linspace(-1, 1, 2*n_elements + 1)';  % 9 nodes
[x_nodes,III]  = sort(xi) ;
n_nodes = length(x_nodes);
n_elements = 1 + (n_nodes-3)/2 ;

% 2. Define connectivity (node indices per element)
elem_nodes = zeros(n_elements, 3);
for e = 1:n_elements
    elem_nodes(e, :) = [2*e-1, 2*e, 2*e+1];
end

% 3. Gauss quadrature (3-point for exact integration of degree 5)
xi_q = [-sqrt(3/5), 0, sqrt(3/5)];
w_q = [5/9, 8/9, 5/9];

% 4. Reference shape functions
phi = {@(xi) 0.5*xi.*(xi - 1), ...
    @(xi) (1 - xi.^2), ...
    @(xi) 0.5*xi.*(xi + 1)};

% 5. Assemble global mass matrix
if ORTHOGONALITY_L2 ==1
    M = zeros(n_nodes);
    for e = 1:n_elements
        nodes = elem_nodes(e, :);
        xe = x_nodes(nodes);
        Me = zeros(3, 3);
        for q = 1:length(xi_q)
            xi = xi_q(q);
            w = w_q(q);
            N = cellfun(@(f) f(xi), phi);
            J = (xe(3) - xe(1)) / 2;
            Me = Me + (N' * N) * w * J;
        end
        % Assemble local Me into global M
        for i = 1:3
            for j = 1:3
                M(nodes(i), nodes(j)) = M(nodes(i), nodes(j)) + Me(i, j);
            end
        end
    end
else
    M = speye(n_nodes,n_nodes) ;
end

% 6. Build monomial basis at global nodes
%max_degree = n_nodes - 1;  % max degree = 8
V = zeros(n_nodes, max_degree + 1);
for k = 0:max_degree
    V(:, k+1) = x_nodes.^k;
end

% 7. Orthogonalize with respect to M (with reorthogonalization)
Q = zeros(size(V));
for i = 1:max_degree+1
    v = V(:, i);
    % First pass
    for j = 1:i-1
        proj = (Q(:, j)' * M * v);
        v = v - proj * Q(:, j);
    end
    % Second pass
    for j = 1:i-1
        proj = (Q(:, j)' * M * v);
        v = v - proj * Q(:, j);
    end
    % Normalize
    normM = sqrt(v' * M * v);
    Q(:, i) = v / normM;
end

% 8. Check orthogonality
G = Q' * M * Q;
% disp('Inner product matrix (should be identity):');
% disp(G);
fprintf('Max off-diagonal error: %.2e\n', max(max(abs(G - eye(size(G))))));



if PLOT_FUNCTIONS == 1
    % 9. Plot orthogonal polynomials
    figure;
    hold on;
    for k = 1:size(Q,2)
        plot(x_nodes, Q(:,k), 'o-', 'DisplayName', sprintf('P_%d', k-1));
    end
    legend show;
    xlabel('x');
    ylabel('P_n(x)');
    title('Orthogonal polynomials on FEM quadratic mesh (4 elements)');
    grid on;
    
    
    Q = Q(III,:) ;
    M = M(III,III);
    
end
