classdef Geometry_Surface < Geometry
    
    properties (Access = public)
        normalVector
    end
    
    properties (Access = private)
       drDtxi 
       jacobian
    end
    
    methods (Access = public)
        
        function obj = Geometry_Surface(cParams)
            obj.init(cParams);
        end

        function detJ = computeJacobianDeterminant(obj,xV)
            J = obj.computeJacobian(xV);
            n = obj.computeNormalVectors(J);
            detJ = squeeze(pagenorm(n));
        end

    end
    
    methods (Access = private)

        function J = computeJacobian(obj,xV)
            nDimGlo  = size(obj.coord,1);
            nElem    = size(obj.coord,3);
            dShapes  = obj.xFE.computeShapeDerivatives(xV);
            nDimElem = size(dShapes,1);
            nPoints  = size(xV,2);
            J = zeros(nDimElem,nDimGlo,nPoints,nElem);
            for iDimGlo = 1:nDimGlo
                for iDimElem = 1:nDimElem
                        dShapeIK = dShapes(iDimElem,:,:);
                        dShapeIK = squeeze(dShapeIK);
                        xKJ      = permute(obj.xFE.fValues(iDimGlo,:,:),[2 3 1]);
                        xKJ      = squeeze(xKJ);
                        jacIJ    = dShapeIK*xKJ;
                        J(iDimElem,iDimGlo,:,:) = squeeze(permute(J(iDimElem,iDimGlo,:,:),[3 4 1 2])) + jacIJ;
                end
            end
        end
        
        function normalVector = computeNormalVectors(obj,J)
            nDimGlo = size(J,2);
            nPoints = size(J,3);
            nElem = size(J,4);

            normalVector = zeros(1,nDimGlo,nPoints,nElem);
            DxDxi  = J(1,1,:,:);
            DxDeta = J(2,1,:,:);
            DyDxi  = J(1,2,:,:);
            DyDeta = J(2,2,:,:);
            DzDxi  = J(1,3,:,:);
            DzDeta = J(2,3,:,:);
            normalVector(:,1,:,:) = DyDxi.*DzDeta - DzDxi.*DyDeta;
            normalVector(:,2,:,:) = DxDxi.*DzDeta - DzDxi.*DxDeta;
            normalVector(:,3,:,:) = DxDxi.*DyDeta - DyDxi.*DxDeta;
        end

    end
end