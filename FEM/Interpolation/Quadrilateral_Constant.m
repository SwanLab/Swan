classdef Quadrilateral_Constant < Interpolation

    methods (Access = public)

        function obj = Quadrilateral_Constant(cParams)
            obj.init(cParams);
            obj.computeParams();
        end

        function computeShapeDeriv(obj,posgp)
            obj.computeShapes()
            obj.computeShapeDerivatives();
        end

    end

    methods (Access = private)

        function computeParams(obj)
            obj.type      = 'QUAD';
            obj.ndime     = 2;
            obj.nnode     = 1;
            obj.pos_nodes = [0 0];
        end

        function computeShapes(obj)
            obj.shape = [1];

        end

        function computeShapeDerivatives(obj)
            obj.deriv = [0; 0];
        end

    end
end