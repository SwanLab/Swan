classdef Quadrilateral_Bilinear < Interpolation
    
    methods (Access = public)
        
        function obj = Quadrilateral_Bilinear(cParams)
            obj.init(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function computeParams(obj)
            obj.type = 'QUADRILATERAL';
            obj.ndime = 2;
            obj.nnode = 4;
            obj.pos_nodes = [-1 -1; 1 -1; 1 1; -1 1];
            % obj.isoDv = 4;
        end

        function shape = evaluateShapeFunctions(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            s = xV(1,:,:);
            t = xV(2,:,:);
            I = ones(size(t));
            shape = zeros(obj.nnode,ngaus,nelem);
            shape(1,:,:) = 0.25*(I-t-s+s.*t);
            shape(2,:,:) = 0.25*(I-t+s-s.*t);
            shape(3,:,:) = 0.25*(I+t+s+s.*t);
            shape(4,:,:) = 0.25*(I+t-s-s.*t);
        end
        
        function deriv = evaluateShapeDerivatives(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            s = xV(1,:,:);
            t = xV(2,:,:);
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

end
