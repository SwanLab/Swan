classdef Geometry_Surface < Geometry
    
    methods (Access = public)
        
        function obj = Geometry_Surface(cParams)
            obj.init(cParams);
        end

        function detJ = computeJacobianDeterminant(obj,xV)
            J = obj.computeJacobian(xV);
            n = obj.computeNormalVectors(J);
            detJ = squeeze(pagenorm(n));
        end

        function n = computeNormals(obj, xV)
            J = obj.computeJacobian(xV);
            n = obj.computeNormalVectors(J);
        end
        

    end
    
    methods (Access = private)
        
        function normalVector = computeNormalVectors(obj,J)
            nDimGlo = size(J,2);
            nPoints = size(J,3);
            nElem = size(J,4);

            normalVector = zeros(1,nDimGlo,nPoints,nElem);
            DxDxi  = squeeze(J(1,1,:,:))';
            DxDeta = squeeze(J(2,1,:,:))';
            DyDxi  = squeeze(J(1,2,:,:))';
            DyDeta = squeeze(J(2,2,:,:))';
            DzDxi  = squeeze(J(1,3,:,:))';
            DzDeta = squeeze(J(2,3,:,:))';
            normalVector(:,1,:,:) = DyDxi.*DzDeta - DzDxi.*DyDeta;
            normalVector(:,2,:,:) = DzDxi.*DxDeta - DxDxi.*DzDeta;
            normalVector(:,3,:,:) = DxDxi.*DyDeta - DyDxi.*DxDeta;
        end

    end
end