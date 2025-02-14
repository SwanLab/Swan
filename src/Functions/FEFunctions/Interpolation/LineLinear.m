classdef LineLinear < Interpolation
    
    methods (Access = public)
        
        function obj = LineLinear(cParams)
            obj.init(cParams);
        end
        
    end
    
    methods (Access = protected)

        function computeParams(obj)
            obj.ndime = 1;
            obj.nnode = 2;
            obj.pos_nodes = [-1; 1];
        end

        function shape = evaluateShapeFunctions(obj, xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            s = xV(1,:,:);
            I = ones(size(s));
            shape = zeros(obj.nnode,ngaus,nelem);
            shape(1,:,:) = 0.5*(I-s);
            shape(2,:,:) = 0.5*(s+I);
        end

        function deriv = evaluateShapeDerivatives(obj, xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            deriv = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            deriv(1,1,:,:) = -0.5;
            deriv(1,2,:,:) = 0.5;
        end
        
    end
    
    methods (Access = public)
        function deriv2 = evaluateShapeSecondDerivatives(obj, xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            deriv2 = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            deriv2(1,1,:,:) = 0;
            deriv2(1,2,:,:) = 0;
        end
    end
end
