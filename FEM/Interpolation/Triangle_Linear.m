classdef Triangle_Linear < Interpolation
    
    methods (Access = public)
        
        function obj = Triangle_Linear(cParams)
            obj.init(cParams);
            obj.computeParams();
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
            % obj.isoDv = 0.5;
        end
        
    end
    
end
