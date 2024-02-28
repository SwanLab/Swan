classdef MatrixVectorizedInverterFactory < handle
    
    methods (Static, Access = public)
        
        function inverter = create(A)
            ndime = size(A,1);
            switch ndime
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
            ndime = size(A,1);
        end
        
    end
    
end