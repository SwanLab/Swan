classdef Line_Quadratic < Interpolation
    
    methods (Access = public)
        
        function obj = Line_Quadratic(cParams)
            obj.init(cParams);
        end
        
    end
    
    methods (Access = protected)

        function computeParams(obj)
            obj.ndime = 1;
            obj.nnode = 3;
            obj.pos_nodes = [-1; 1; 0];
        end
        
        function shape = evaluateShapeFunctions(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            s = xV(1,:,:);
            shape = zeros(obj.nnode,ngaus,nelem);
            shape(1,:,:) = (s.*(s - 1))./2;
            shape(2,:,:) = (s.*(s + 1))./2;
            shape(3,:,:) = -(s - 1).*(s + 1);
        end
        
        function deriv = evaluateShapeDerivatives(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            s = xV(1,:,:);
            deriv = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            deriv(1,1,:,:) = s-0.5;
            deriv(1,2,:,:) = s+0.5;
            deriv(1,3,:,:) = -2.*s;
        end
        
    end
    
end
