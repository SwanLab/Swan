classdef RebuildKMTData < handle
    % process K, M, T data tu use in EIFEM
    
    properties (Access = private)
        
    end
    
    methods (Static, Access = public)
        
        function [K,M,T] = compute(method,r,mesh)

            % Initialize
            K = []; M = []; T = [];

            %K: Kcoarse Prediction from NN after Chol
            %M: Mcoarse Prediction from NN after Chol
            %T: T prediction from NN
            load('Kcoarse_predictor.mat');
            load('McoarsePredictorNN.mat');

            Kc = opt.computeOutputValues([r]);
            Mc = McPredictorNN.computeOutputValues([r]);

            K = RebuildKMTData.rebuildMatrix(Kc);
            M = RebuildKMTData.rebuildMatrix(Mc);

            if(method=='A')
                load('TpredictorNN.mat');
                T = RebuildKMTData.rebuildTPrediction(r,mesh,tPredictorNN);
            % elseif(method=='B')
            %     T = rebuildTMethodB(U,S,V);
            % else
            %     T = rebuildTMethodC(U,S,V);
            end

        end
        
    end
    
    methods (Static, Access = private)

        function KMmatrix = rebuildMatrix(KM_NN)
            % Processes K or M prediction from NN and outputs a 8x8
            % Kcoarse or Mcoarse matrix.
            
            n = 8; %matrix dimensinos
            L = zeros(n);
            idx = tril(true(n));
            L(idx) = KM_NN; % fill lower triangle
            d = diag(L);
            d(d <= 0) = eps;                      
            L(1:n+1:end) = d;
            KMmatrix = L*L.';%revert cholesky decomposition   
            
        end
        
        function T = rebuildTPrediction(r, mesh, NN)
            
            % Initialize empty matrix (or pre-allocate for speed if possible)
            Taux2 = [];
            
            for i = 1:size(mesh.coord, 1)  % Evaluates all coordinates
                dataInput = [r, mesh.coord(i,:)];  
                
                % Format: [Tx1, Ty1, Tx2, Ty2, ..., Tx8, Ty8]
                T_NN = NN.computeOutputValues(dataInput);
                
                % Reshape into 2x8
                T_block = reshape(T_NN, 2, 8);
                
                % 3. Stack vertically
                % Final T will have dimensions [2*N_coords x 8]
                Taux2 = cat(1, Taux2, T_block);
            end
            
            T = Taux2;
        
        end
                
    end
    
end
