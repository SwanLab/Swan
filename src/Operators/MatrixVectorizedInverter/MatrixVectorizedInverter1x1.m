classdef MatrixVectorizedInverter1x1 < MatrixVectorizedInverterInterface
    
    methods (Access = public)
        
        function B = computeInverse(~,A)
            B = 1./A;
        end
        
        function det = computeDeterminant(~,A)
            det = A;
        end
        
    end
    
end