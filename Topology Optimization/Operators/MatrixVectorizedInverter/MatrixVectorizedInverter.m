classdef MatrixVectorizedInverter < MatrixVectorizedInverter_Interface
    
    properties (Access = private)
        executor
    end
    
    methods (Access = public)
        
        function inv = computeInverse(obj,A)
            obj.createExecutor(A);
            inv = obj.executor.computeInverse(A);
        end
        
        function det = computeDeterminant(obj,A)
            obj.createExecutor(A);
            det = obj.executor.computeDeterminant(A);
        end
        
    end
    
    methods (Access = private)
        
        function createExecutor(obj,A)
            obj.executor = MatrixVectorizedInverterFactory().create(A);
        end
        
    end
    
end