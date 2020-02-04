classdef Triangle_Linear < Interpolation
    
    methods (Access = public)
        
        function obj = Triangle_Linear(cParams)
            obj.init(cParams);
            obj.computeParams();
            obj.computeCases();
            obj.computeCoordAndConnec();
        end
        
        function computeShapeDeriv(obj,posgp)
            obj.computeShapes(posgp);
            obj.computeShapeDerivatives(posgp);
        end
        
    end
    
    methods (Access = private)
        
        function computeParams(obj)
            obj.type  = obj.mesh.geometryType;
            obj.ndime = 2;
            obj.nnode = 3;
            obj.pos_nodes = [0 0; 1 0; 0 1];
            obj.isoDv = 0.5;
            obj.main_loop = [3 3];
            obj.extra_cases = [];                        
        end
        
        function computeShapes(obj,posgp)
            ngaus = size(posgp,2);
            s = posgp(1,:);
            t = posgp(2,:);
            I = ones(1,ngaus);            
            obj.shape = [I-s-t;s;t];
        end
        
        function computeShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);                      
            derivative = [-1.0 1.0 0.0;-1.0 0.0 1.0];                        
            obj.deriv = repmat(derivative,1,1,ngaus);            
        end
        
        function computeCases(obj)
            obj.cases(:,:,1) = [1 4 5;
                4 2 3;
                5 4 3];
            obj.cases(:,:,2) = [1 4 3;
                4 2 5;
                4 5 3];
            obj.cases(:,:,3) = [1 4 5;
                1 2 4;
                5 4 3];
            obj.selectcases =  [1     0;
                2     0;
                3     3;
                0     2;
                0     1];
        end
        
    end
    
end
