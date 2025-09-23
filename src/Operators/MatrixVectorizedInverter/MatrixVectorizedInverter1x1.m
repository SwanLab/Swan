classdef MatrixVectorizedInverter1x1 < MatrixVectorizedInverterInterface
    
    methods (Access = public)
        
        function B = computeInverse(~,A)
            B = 1./A;
        end
        
        function det = computeDeterminant(~,A)
            det = squeeze(pagenorm(A));
        end
        
    end
    
end