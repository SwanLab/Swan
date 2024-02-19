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

    methods (Access = private)

        function J = computeJacobian(obj,xV)
            nDimGlo  = size(obj.coord,1);
            nElem    = size(obj.coord,3);
            dShapes  = obj.interpolation.computeShapeDerivatives(xV);
            nDimElem = size(dShapes,1);
            nPoints  = size(xV,2);
            J = zeros(nDimElem,nDimGlo,nPoints,nElem);
            for iDimGlo = 1:nDimGlo
                for iDimElem = 1:nDimElem
                        dShapeIK = squeezeParticular(dShapes(iDimElem,:,:),1)';
                        xKJ = squeezeParticular(obj.coord(iDimGlo,:,:),1);
                        jacIJ    = dShapeIK*xKJ;
                        J(iDimElem,iDimGlo,:,:) = squeezeParticular(J(iDimElem,iDimGlo,:,:),[1 2]) + jacIJ;
                end
            end
        end

    end

end