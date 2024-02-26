classdef Tetrahedra_Linear < Interpolation
    
    methods (Access = public)
        
        function obj = Tetrahedra_Linear(cParams)
            obj.init(cParams);
            obj.computeParams();
        end

        function shape = computeShapeFunctions(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            s = posgp(1,:,:);
            t = posgp(2,:,:);
            u = posgp(3,:,:);
            shape = zeros(obj.nnode,ngaus,nelem);
            shape(1,:,:) = 1-t-s-u;
            shape(2,:,:) = s;
            shape(3,:,:) = t;
            shape(4,:,:) = u;
        end

        function deriv = computeShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            deriv = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            deriv(1,1,:,:) = -1;
            deriv(1,2,:,:) = 1;
            deriv(1,3,:,:) = 0;
            deriv(1,4,:,:) = 0;
            deriv(2,1,:,:) = -1;
            deriv(2,2,:,:) = 0;
            deriv(2,3,:,:) = 1;
            deriv(2,4,:,:) = 0;
            deriv(3,1,:,:) = -1;
            deriv(3,2,:,:) = 0;
            deriv(3,3,:,:) = 0;
            deriv(3,4,:,:) = 1;
        end
        
    end
    
    methods (Access = private)
        
        function computeParams(obj)
            obj.ndime = 3;
            obj.nnode = 4;
            obj.pos_nodes = [0 0 0;
                1 0 0;
                0 1 0;
                0 0 1];
            % obj.isoDv = 1/6;
        end
    end
    
end
