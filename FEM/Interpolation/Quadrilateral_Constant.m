classdef Quadrilateral_Constant < Interpolation

    methods (Access = public)

        function obj = Quadrilateral_Constant(cParams)
            obj.init(cParams);
            obj.computeParams();
        end
        
        function shape = computeShapeFunctions(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            N = ones(obj.nnode,ngaus,nelem);
            shape = N;
        end
        
        function deriv = computeShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            dN = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            deriv = dN;
        end

    end

    methods (Access = private)

        function computeParams(obj)
            obj.type      = 'QUAD';
            obj.ndime     = 2;
            obj.nnode     = 1;
            obj.pos_nodes = [0 0];
        end

    end
end