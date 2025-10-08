classdef StiffnessMatrixHandler < handle
    
    methods (Static)
        
        function K = vectorToMatrix(vector)
            % Convierte un vector de 64 elementos a una matriz 8x8
            % INPUT: vector - vector de 64 elementos
            % OUTPUT: K - matriz 8x8 de rigidez
            K = reshape(vector, 8, 8);
        end
        
        function vector = matrixToVector(K)
            % Convierte una matriz 8x8 a un vector de 64 elementos
            % INPUT: K - matriz 8x8 de rigidez
            % OUTPUT: vector - vector de 64 elementos
            vector = reshape(K, 1, 64);
        end
        
        function isValid = validateStiffnessMatrix(K, tolerance)
            % Valida que una matriz sea una matriz de rigidez válida
            % INPUT: K - matriz 8x8
            %        tolerance - tolerancia para comparaciones numéricas
            % OUTPUT: isValid - booleano indicando si es válida
            
            if nargin < 2
                tolerance = 1e-10;
            end
            
            isValid = true;
            
            % Verificar que sea cuadrada
            if size(K, 1) ~= 8 || size(K, 2) ~= 8
                isValid = false;
                return;
            end
            
            % Verificar simetría
            if ~issymmetric(K, tolerance)
                isValid = false;
                return;
            end
            
            % Verificar que sea definida positiva
            try
                eigenvals = eig(K);
                if any(eigenvals <= tolerance)
                    isValid = false;
                    return;
                end
            catch
                isValid = false;
                return;
            end
        end
        
        function K_symmetric = enforceSymmetry(K)
            % Fuerza la simetría en una matriz
            % INPUT: K - matriz 8x8
            % OUTPUT: K_symmetric - matriz simétrica
            K_symmetric = (K + K') / 2;
        end
        
        function K_positive = enforcePositiveDefinite(K, min_eigenval)
            % Fuerza que una matriz sea definida positiva
            % INPUT: K - matriz 8x8
            %        min_eigenval - valor mínimo para los eigenvalores
            % OUTPUT: K_positive - matriz definida positiva
            
            if nargin < 2
                min_eigenval = 1e-6;
            end
            
            % Hacer simétrica primero
            K_sym = StiffnessMatrixHandler.enforceSymmetry(K);
            
            % Calcular eigenvalores y eigenvectores
            [V, D] = eig(K_sym);
            
            % Asegurar que todos los eigenvalores sean positivos
            eigenvals = diag(D);
            eigenvals(eigenvals < min_eigenval) = min_eigenval;
            
            % Reconstruir la matriz
            K_positive = V * diag(eigenvals) * V';
        end
        
        function K_valid = makeValidStiffnessMatrix(K)
            % Convierte una matriz en una matriz de rigidez válida
            % INPUT: K - matriz 8x8
            % OUTPUT: K_valid - matriz de rigidez válida (simétrica y definida positiva)
            
            % Enforzar simetría
            K_sym = StiffnessMatrixHandler.enforceSymmetry(K);
            
            % Enforzar definida positiva
            K_valid = StiffnessMatrixHandler.enforcePositiveDefinite(K_sym);
        end
        
        function plotStiffnessMatrix(K, title_str)
            % Visualiza una matriz de rigidez
            % INPUT: K - matriz 8x8
            %        title_str - título del gráfico
            
            if nargin < 2
                title_str = 'Stiffness Matrix';
            end
            
            figure;
            imagesc(K);
            colorbar;
            title(title_str);
            xlabel('DOF');
            ylabel('DOF');
            
            % Agregar valores en cada celda
            for i = 1:8
                for j = 1:8
                    text(j, i, sprintf('%.2e', K(i,j)), ...
                         'HorizontalAlignment', 'center', ...
                         'FontSize', 8);
                end
            end
        end
        
        function [condition_number, max_eigenval, min_eigenval] = analyzeMatrix(K)
            % Analiza las propiedades de una matriz de rigidez
            % INPUT: K - matriz 8x8
            % OUTPUT: condition_number - número de condición
            %         max_eigenval - máximo eigenvalor
            %         min_eigenval - mínimo eigenvalor
            
            eigenvals = eig(K);
            max_eigenval = max(eigenvals);
            min_eigenval = min(eigenvals);
            condition_number = max_eigenval / min_eigenval;
        end
        
        function cost = stiffnessMatrixCost(K, target_properties)
            % Calcula el costo de una matriz de rigidez basado en propiedades objetivo
            % INPUT: K - matriz 8x8
            %        target_properties - estructura con propiedades objetivo
            % OUTPUT: cost - valor del costo
            
            if nargin < 2
                target_properties = struct();
            end
            
            cost = 0;
            
            % Costo por no ser simétrica
            if isfield(target_properties, 'enforce_symmetry') && target_properties.enforce_symmetry
                symmetry_error = norm(K - K', 'fro');
                cost = cost + target_properties.symmetry_weight * symmetry_error;
            end
            
            % Costo por no ser definida positiva
            if isfield(target_properties, 'enforce_positive_definite') && target_properties.enforce_positive_definite
                eigenvals = eig(K);
                negative_eigenvals = eigenvals(eigenvals < 0);
                if ~isempty(negative_eigenvals)
                    cost = cost + target_properties.positive_definite_weight * sum(negative_eigenvals.^2);
                end
            end
            
            % Costo por desviación de propiedades específicas
            if isfield(target_properties, 'target_diagonal')
                diagonal_error = norm(diag(K) - target_properties.target_diagonal);
                cost = cost + target_properties.diagonal_weight * diagonal_error;
            end
        end
        
    end
end
