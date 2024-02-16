classdef Geometry_Volumetric < Geometry
    
    properties (GetAccess = public, SetAccess = private)
        dNdx
    end
    
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
            coord = obj.xFE.fValues;
            nDime   = size(coord,1);
            nNode   = size(coord,2);
            nElem   = size(coord,3);
            nPoints = size(xV,2);
            dShapes = obj.xFE.computeShapeDerivatives(xV);
            J = zeros(nDime,nDime,nPoints,nElem);
            for iPoints = 1:nPoints
                for iNode = 1:nNode
                    dShapeIK = repmat(dShapes(:,iNode,iPoints),[1 1 nElem]);
                    xKJ      = pagetranspose(obj.xFE.fValues(:,iNode,:));
                    jacIJ    = pagemtimes(dShapeIK, xKJ);
                    J(:,:,iPoints,:) = squeeze(J(:,:,iPoints,:)) + jacIJ;
                end
            end
        end
        
    end
    
end