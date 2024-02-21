classdef Geometry_Volumetric < Geometry

    properties (Access = private)
        matrixInverter
    end

    methods (Access = public)

        function obj = Geometry_Volumetric(cParams)
            obj.init(cParams);
            obj.matrixInverter = MatrixVectorizedInverter();
        end

        function detJ = computeJacobianDeterminant(obj,xV)
            J = obj.computeJacobian(xV);
            detJ = obj.matrixInverter.computeDeterminant(J);
        end

        function invJ = computeInverseJacobian(obj,xV)
            J = obj.computeJacobian(xV);
            invJ = obj.matrixInverter.computeInverse(J);
        end

    end

end