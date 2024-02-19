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
        
    end
    
    
end