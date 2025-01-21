classdef MatrixVectorizedInverter_1x1 < MatrixVectorizedInverter_Interface
    
    methods (Access = public)
        
        function B = computeInverse(~,A)
            B = 1./A;
        end
        
        function det = computeDeterminant(~,A)
            det = A;
        end
        
    end
    
end