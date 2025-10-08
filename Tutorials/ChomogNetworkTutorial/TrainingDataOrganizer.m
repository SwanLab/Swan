classdef TrainingDataOrganizer
    % SPDDataProcessor  Processes stiffness SPD matrices for NN training
    %   Takes a cell array of 8x8 SPD matrices and their labels,
    %   applies Cholesky decomposition, and outputs a numeric matrix
    %   suitable for training neural networks.
    
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
                
                % Check SPD via Cholesky
                [L, p] = chol(A, 'lower');
                if p > 0
                    error('Matrix %d is not SPD.', i);
                end
                
                % Extract lower-triangular part as feature vector
                features = L(tril(true(size(L))))';
                
                % Combine with label
                trainingMatrix(i, :) = [obj.Labels(i), features];
            end
        end
    end
end
