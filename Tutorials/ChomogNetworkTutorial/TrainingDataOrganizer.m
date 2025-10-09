classdef TrainingDataOrganizer
    % TrainingDataOrganizer  Processes stiffness SPD matrices for NN training
    
    properties
        Matrices   % Cell array of SPD matrices
        Labels     % Numeric vector of labels
    end
    
    methods

        function obj = TrainingDataOrganizer(labels, matrices)
            if nargin > 0
                if ~iscell(matrices)
                    error('Object with K matrices must be a cell array of matrices.');
                end
                if numel(matrices) ~= numel(labels)
                    error('Number of matrices and labels must match.');
                end
                obj.Matrices = matrices;
                obj.Labels = labels;
            end
        end
        
        % Generate training matrix
        function trainingMatrix = generateTrainingMatrix(obj)
            n = numel(obj.Matrices);
            
            % Example Cholesky of first matrix to get size
            L = chol(obj.Matrices{1}, 'lower');
            featureLength = numel(L(tril(true(size(L))))); % lower-triangular
            
            trainingMatrix = zeros(n, featureLength + 1);
            
            for i = 1:n
                A = obj.Matrices{i};
                
                                [L, p] = chol(A, 'lower'); % Check SPD via Cholesky decomsd
                if p > 0
                    error('Matrix %d is not SPD.', i);
                end
                
                % Extract lower-triangular part
                features = L(tril(true(size(L))))';
                                
                trainingMatrix(i, :) = [obj.Labels(i), features]; % Combine with label
            end
        end
    end
end
