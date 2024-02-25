classdef Quadrilateral_Bilinear < Interpolation
    
    methods (Access = public)
        
        function obj = Quadrilateral_Bilinear(cParams)
            obj.init(cParams);
            obj.computeParameters();
        end
        function shape = computeShapeFunctions(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            s = posgp(1,:,:);
            t = posgp(2,:,:);
            I = ones(size(t));
            shape = zeros(obj.nnode,ngaus,nelem);
            shape(1,:,:) = 0.25*(I-t-s+s.*t);
            shape(2,:,:) = 0.25*(I-t+s-s.*t);
            shape(3,:,:) = 0.25*(I+t+s+s.*t);
            shape(4,:,:) = 0.25*(I+t-s-s.*t);
        end
        
        function deriv = computeShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            s = posgp(1,:,:);
            t = posgp(2,:,:);
            I = ones(size(t));
            deriv = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            deriv(1,1,:,:) = 0.25*(-I+t);
            deriv(1,2,:,:) = 0.25*(+I-t);
            deriv(1,3,:,:) = 0.25*(+I+t);
            deriv(1,4,:,:) = 0.25*(-I-t);
            deriv(2,1,:,:) = 0.25*(-I+s);
            deriv(2,2,:,:) = 0.25*(-I-s);
            deriv(2,3,:,:) = 0.25*(+I+s);
            deriv(2,4,:,:) = 0.25*(+I-s);
        end
        
    end
    
    methods (Access = private)
        
        function computeParameters(obj)
            obj.type = 'QUADRILATERAL';
            obj.ndime = 2;
            obj.nnode = 4;
            obj.pos_nodes = [-1 -1; 1 -1; 1 1; -1 1];
            % obj.isoDv = 4;
        end
        
    end

end
