classdef Tetrahedra_Constant < Interpolation

    properties (Access = private)
        ngaus
    end
       
    methods (Access = public)

        function obj = Tetrahedra_Constant(cParams)
            obj.init(cParams);
            obj.computeParams();
        end

        function computeShapeDeriv(obj,posgp)
            obj.ngaus = size(posgp,2);
            obj.computeShapes()
            obj.computeShapeDerivatives();
        end

    end

    methods (Access = private)

        function computeParams(obj)
            obj.type = 'TETRAHEDRA';
            obj.ndime = 3;
            obj.nnode = 1;
            obj.pos_nodes = [1/4 1/4 1/4];
        end

        function computeShapes(obj)
            N = ones(obj.nnode,obj.ngaus);
            obj.shape = N;

        end

        function computeShapeDerivatives(obj)
            dN = zeros(obj.ndime,obj.nnode,obj.ngaus);
            obj.deriv = dN;
        end

    end
end
