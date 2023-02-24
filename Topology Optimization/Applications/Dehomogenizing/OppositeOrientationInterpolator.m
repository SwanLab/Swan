classdef OppositeOrientationInterpolator < handle
    
    
    properties (Access = private)
       mesh
       orientationVector
    end
    
    methods (Access = public)
        
        function obj = OppositeOrientationInterpolator(cParams)
            obj.init(cParams);
        end
        
        function c = compute(obj)
            c = obj.computeSymmetricCondition();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh              = cParams.mesh;
            obj.orientationVector = cParams.orientationVector;
        end
        
        function sC = computeSymmetricCondition(obj)
            nnodeD    = obj.mesh.nnodeElem;
            nElemD    = obj.mesh.nelem;
            nnodesC   = obj.mesh.nnodes;
            connecC   = obj.mesh.connec;
            connecD   = obj.computeDiscontinuousConnectivities();
            isCoherent = obj.orientationVector.isCoherent;
            sC = sparse(nnodeD*nElemD,nnodesC);
            for iNode = 1:nnodeD
                isC  = squeeze(isCoherent.fValues(1,iNode,:));
                cond = obj.computeConformalMapCondition(isC);
                nodesC = connecC(:,iNode);
                nodesD = connecD(:,iNode);
                sC = sC + sparse(nodesD,nodesC,cond,nnodeD*nElemD,nnodesC);
            end
        end
       
        function connec = computeDiscontinuousConnectivities(obj)
            nNodes = obj.mesh.nnodeElem*obj.mesh.nelem;
            nodes  = 1:nNodes;
            connec = reshape(nodes,obj.mesh.nnodeElem,obj.mesh.nelem)';
        end
    end

    methods (Access = private, Static)

        function m = computeConformalMapCondition(isCoherent)
            isMapSymmetric = isCoherent;
            isMapAntiSymm  = ~isCoherent;
            m = zeros(size(isMapSymmetric));
            m(isMapSymmetric) = 1;
            m(isMapAntiSymm)  = -1;
        end

    end
    
end