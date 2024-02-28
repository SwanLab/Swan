classdef Line_Quadratic < Interpolation
    
    methods (Access = public)
        
        function obj = Line_Quadratic(cParams)
            obj.init(cParams);
            obj.computeParams();
            obj.computeCases();
        end
        
        function computeShapeDeriv(obj,posgp)
            obj.computeShapes(posgp);
            obj.computeShapeDerivatives(posgp);
        end
        
    end
    
    methods (Access = private)

        function computeParams(obj)
            obj.ndime = 1;
            obj.nnode = 3;
            obj.pos_nodes = [-1; 1; 0];
            obj.isoDv = 2;
        end
        
        function computeShapes(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            s = posgp(1,:,:);
            obj.shape = zeros(obj.nnode,ngaus,nelem);
            obj.shape(1,:,:) = (s.*(s - 1))./2;
            obj.shape(2,:,:) = (s.*(s + 1))./2;
            obj.shape(3,:,:) = -(s - 1).*(s + 1);
        end
        
        function computeShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            s = posgp(1,:,:);
            obj.deriv = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            obj.deriv(1,1,:,:) = s-0.5;
            obj.deriv(1,2,:,:) = s+0.5;
            obj.deriv(1,3,:,:) = -2.*s;
        end
        
        function computeCases(obj)
            obj.iteration = [1; 
                             2 ];
        end
        
    end
    
end
