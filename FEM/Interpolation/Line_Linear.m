classdef Line_Linear < Interpolation
    
    methods (Access = public)
        
        function obj = Line_Linear(cParams)
            obj.init(cParams);
            obj.computeParams();
        end
        
        function shape = computeShapeFunctions(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            s = posgp(1,:,:);
            I = ones(size(s));
            shape = zeros(obj.nnode,ngaus,nelem);
            shape(1,:,:) = 0.5*(I-s);
            shape(2,:,:) = 0.5*(s+I);
        end
        
        function deriv = computeShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            deriv = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            deriv(1,1,:,:) = -0.5;
            deriv(1,2,:,:) = 0.5;
        end
        
    end
    
    methods (Access = private)

        function computeParams(obj)
            obj.ndime = 1;
            obj.nnode = 2;
            obj.pos_nodes = [-1; 1];
            % obj.isoDv = 2;
        end
        
    end
    
end
