classdef MatrixVectorizedInverter < MatrixVectorizedInverterInterface
    
    properties (Access = private)
        executor
    end
    
    methods (Static, Access = public)
        
        function inv = computeInverse(A)
            exec = MatrixVectorizedInverterFactory.create(A);
            inv = exec.computeInverse(A);
        end
        
        function det = computeDeterminant(A)
            exec = MatrixVectorizedInverterFactory.create(A);
            det = exec.computeDeterminant(A);
        end
        
    end
    
    methods (Access = private)
        
        function createExecutor(obj,A)
            obj.executor = MatrixVectorizedInverterFactory.create(A);
        end
        
    end
    
end