classdef Line_Cubic < Interpolation
    
    methods (Access = public)
        
        function obj = Line_Cubic(cParams)
            obj.init(cParams);
        end
        
    end
    
    methods (Access = protected)

        function computeParams(obj)
            obj.ndime = 1;
            obj.nnode = 4;
            obj.pos_nodes = [-1; 1; -1/3; 1/3];
        end
        
        function shape = evaluateShapeFunctions(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            s = xV(1,:,:);
            shape = zeros(obj.nnode,ngaus,nelem);
            shape(1,:,:) = -(9.*(s - 1).*(s - 1/3).*(s + 1/3))./16;
            shape(2,:,:) = (9.*(s + 1).*(s - 1/3).*(s + 1/3))./16;
            shape(3,:,:) = (27.*(s - 1).*(s + 1).*(s - 1/3))./16;
            shape(4,:,:) = -(27.*(s - 1).*(s + 1).*(s + 1/3))./16;
        end
        
        function deriv = evaluateShapeDerivatives(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            s = xV(1,:,:);
            obj.deriv = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            deriv(1,1,:,:) = (9.*s)/8 - (27.*s.^2)/16 + 1/16;
            deriv(1,2,:,:) = (27.*s.^2)/16 + (9.*s)/8 - 1/16;
            deriv(1,3,:,:) = (81.*s.^2)/16 - (9.*s)/8 - 27/16;
            deriv(1,4,:,:) = 27/16 - (81.*s.^2)/16 - (9.*s)/8;
        end
        
    end
    
end
