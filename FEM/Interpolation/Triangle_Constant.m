classdef Triangle_Constant < Interpolation

    methods (Access = public)
        
        function obj = Triangle_Constant(cParams)
            obj.init(cParams);            
            obj.computeParams();
            obj.computeCoordAndConnec();                        
        end
        
        function computeShapeDeriv(obj,posgp)
            obj.computeShapes()
            obj.computeShapeDerivatives();
        end
        
    end
    
    methods (Access = private)

        function computeParams(obj)
            obj.type = 'TRIANGLE';
            obj.ndime = 2;
            obj.nnode = 1;
            obj.pos_nodes = [1/3 1/3];
        end
        
        function computeShapes(obj)
            obj.shape = @(s,t) {1};            
        end
        
        function computeShapeDerivatives(obj)
            obj.deriv = @(s,t) {0;0};            
        end
 
    end
end
