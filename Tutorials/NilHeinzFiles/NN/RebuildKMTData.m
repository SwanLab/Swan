classdef RebuildKMTData < handle
    % process K, M, T data tu use in EIFEM
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function [K,M,T] = RebuildKMTData(K,M,method,T,U,S,V)
            %K: Kcoarse Prediction from NN after Chol
            %M: Mcoarse Prediction from NN after Chol
            %T: T prediction from NN

            K = rebuildMatrix(K);
            M = rebuildMatrix(M);

            if(method=='A')
                T = rebuildTPrediction(T);
            elseif(method=='B')
                T = rebuildTMethodB(U,S,V);
            else
                T = rebuildTMethodC(U,S,V);
            end

        end
        
    end
    
    methods (Access = private)

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
        
        function T = rebuildTPrediction(T)
            
        end
                
    end
    
end
