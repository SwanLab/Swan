classdef LineCubic < Interpolation
    
    methods (Access = public)
        
        function obj = LineCubic(cParams)
            obj.init(cParams);
            obj.computeParams();
        end
        
        function shape = computeShapeFunctions(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            s = posgp(1,:,:);
            shape = zeros(obj.nnode,ngaus,nelem);
            shape(1,:,:) = -(9.*(s - 1).*(s - 1/3).*(s + 1/3))./16;
            shape(2,:,:) = (9.*(s + 1).*(s - 1/3).*(s + 1/3))./16;
            shape(3,:,:) = (27.*(s - 1).*(s + 1).*(s - 1/3))./16;
            shape(4,:,:) = -(27.*(s - 1).*(s + 1).*(s + 1/3))./16;
        end
        
        function deriv = computeShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            s = posgp(1,:,:);
            deriv = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            deriv(1,1,:,:) = (9.*s)/8 - (27.*s.^2)/16 + 1/16;
            deriv(1,2,:,:) = (27.*s.^2)/16 + (9.*s)/8 - 1/16;
            deriv(1,3,:,:) = (81.*s.^2)/16 - (9.*s)/8 - 27/16;
            deriv(1,4,:,:) = 27/16 - (81.*s.^2)/16 - (9.*s)/8;
        end
        
    end
    
    methods (Access = private)

        function computeParams(obj)
            obj.ndime = 1;
            obj.nnode = 4;
            obj.pos_nodes = [-1; 1; -1/3; 1/3];
            %obj.isoDv = 2;
        end
        
    end
    
end
