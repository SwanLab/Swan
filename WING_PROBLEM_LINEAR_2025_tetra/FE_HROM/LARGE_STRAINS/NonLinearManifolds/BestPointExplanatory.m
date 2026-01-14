function total_error = BestPointExplanatory(A)
% function total_error = BestPointExplanatory(A)
%
% PURPOSE:
%   Computes the "explanatory power" of each column of matrix A by 
%   measuring how well it can represent the remaining columns via 
%   one-dimensional orthogonal projections.
%
%   The function returns a vector where each entry corresponds to the total 
%   squared residual error incurred when projecting all other columns onto 
%   the current one. The column with the lowest total error is considered 
%   the most representative or explanatory.
%
%   This idea is analogous to the greedy selection criterion in the 
%   Reduced Basis Method (RBM), where the most "informative" snapshot is 
%   chosen at each step.
%
% INPUT:
%   - A : [n × m] matrix with linearly independent columns. Typically,
%         columns represent integrand snapshots or modes at selected 
%         quadrature points.
%
% OUTPUT:
%   - total_error : [1 × m] vector where total_error(j) is the accumulated
%                   squared projection error when column j is used to
%                   approximate all others.
%
% AUTHOR:
%   Joaquín A. Hernández, UPC/CIMNE, 2025
 
[n, m] = size(A);
total_error = zeros(1, m);  % Stores total residual error for each candidate column

for j = 1:m
    v = A(:, j);                      % Candidate basis column
    v_norm_sq = v' * v;              % Precompute ||v||^2
    error_j = 0;

    for i = 1:m
        if i ~= j
            u = A(:, i);             % Other column
            coeff = (v' * u) / v_norm_sq;
            proj = coeff * v;        % Projection of u onto v
            residual = u - proj;
            error_j = error_j + norm(residual)^2;  % Squared residual
        end
    end

    total_error(j) = error_j;
end

% Identify best representative column
% [aaaa,bbb] = sort(total_error);
% 
% fprintf('Best representative column is column #%d\n', best_col_index);
