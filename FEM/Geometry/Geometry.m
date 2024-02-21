classdef Geometry < handle

    properties (SetAccess = private, GetAccess = protected)
        coord
        interpolation
    end

    methods (Access = public, Static)

        function obj = create(cParams)
            f = GeometryFactory();
            obj = f.create(cParams);
        end

    end

    methods (Access = protected)

        function init(obj,cParams)
            obj.coord = cParams.coord;
            obj.interpolation = cParams.interp;
        end

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