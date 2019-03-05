classdef MatrixVectorizedInverterFactory < handle
    
    methods (Access = public)
        
        function inverter = create(obj,A)
            switch obj.computeDimension(A)
                case 1
                    inverter = MatrixVectorizedInverter_1x1();
                case 2
                    inverter = MatrixVectorizedInverter_2x2();
                case 3
                    inverter = MatrixVectorizedInverter_3x3();
            end
        end
        
    end
    
    methods (Access = private, Static)
        
        function ndime = computeDimension(A)
            ndime = min(size(A));
        end
        
    end
    
end