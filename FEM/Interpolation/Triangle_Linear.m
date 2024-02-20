classdef Triangle_Linear < Interpolation
    
    methods (Access = public)
        
        function obj = Triangle_Linear(cParams)
            obj.init(cParams);
            obj.computeParams();
            obj.computeCases();
        end
        
        function shape = computeShapeFunctions(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            s = posgp(1,:,:);
            t = posgp(2,:,:);
            I = ones(size(t));
            shape = zeros(obj.nnode,ngaus,nelem);
            shape(1,:,:) = I-s-t;
            shape(2,:,:) = s;
            shape(3,:,:) = t;
        end
        
        function deriv = computeShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            deriv = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            deriv(1,1,:,:) = -1;
            deriv(1,2,:,:) = 1;
            deriv(1,3,:,:) = 0;
            deriv(2,1,:,:) = -1;
            deriv(2,2,:,:) = 0;
            deriv(2,3,:,:) = 1;
        end
        
    end
    
    methods (Access = private)
        
        function computeParams(obj)
            obj.ndime = 2;
            obj.nnode = 3;
            obj.pos_nodes = [0 0; 1 0; 0 1];
            obj.isoDv = 0.5;
            obj.main_loop = [3 3];
            obj.extra_cases = [];
        end
        
        function computeCases(obj)
            obj.iteration = [1 2 3;
                             2 3 1];
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
