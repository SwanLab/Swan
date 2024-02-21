classdef Geometry_Line < Geometry
    
    methods (Access = public)
        
        function obj = Geometry_Line(cParams)
            obj.init(cParams);
        end

        function detJ = computeJacobianDeterminant(obj,xV)
            J = obj.computeJacobian(xV);
            detJ = squeeze(pagenorm(J));
        end
        
        function invJ = computeInverseJacobian(obj,xV)
            detJ = obj.computeJacobianDeterminant(xV);
            invJ = 1./detJ;
        end

    end
    
    
end