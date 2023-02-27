classdef OrientationVectors < handle
    
    properties (GetAccess = public, SetAccess = private)
        value
        isCoherent
        interpolator
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
      mesh
      theta
    end
    
    methods (Access = public)
        
        function obj = OrientationVectors(cParams)
            obj.init(cParams)
            obj.createOrientationVector();
            obj.computeIsOrientationCoherent();
            obj.createInterpolator();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.theta = cParams.theta;
        end
        
        function createOrientationVector(obj)
            a1(:,1) = cos(obj.theta);
            a1(:,2) = sin(obj.theta);
            a2(:,1) = -sin(obj.theta);
            a2(:,2) = cos(obj.theta);
            a(:,:,1) = a1;
            a(:,:,2) = a2;
            nDim = obj.mesh.ndim;
            orientation = cell(nDim,1);
            for iDim = 1:nDim
                s.fValues = a(:,:,iDim);
                s.mesh   = obj.mesh;
                bf = P1Function(s);
                orientation{iDim} = bf;
            end 
            obj.value = orientation;
        end
        
        function computeIsOrientationCoherent(obj)
            nnode = obj.mesh.nnodeElem;
            nElem = obj.mesh.nelem;
            isCoh = false(1,nnode,nElem);
            a1D   = obj.value{1}.project('P1D'); 
            a1    = a1D.fValues;
            aN1   = squeeze(a1(:,1,:));
            for iNode = 1:nnode
                aNi    = squeeze(a1(:,iNode,:));
                aN1aNI = dot(aN1,aNi);
                isCoh(1,iNode,:) = (aN1aNI)>0;
            end
            s.fValues = isCoh;
            s.mesh    = obj.mesh;
            isCF       = P1DiscontinuousFunction(s);            
            obj.isCoherent = isCF;
        end

        function sC = createInterpolator(obj)
            nnodeD    = obj.mesh.nnodeElem;
            nElemD    = obj.mesh.nelem;
            nnodesC   = obj.mesh.nnodes;
            connecC   = obj.mesh.connec;
            connecD   = obj.computeDiscontinuousConnectivities();
            isCo = obj.isCoherent;
            sC = sparse(nnodeD*nElemD,nnodesC);
            for iNode = 1:nnodeD
                isC  = squeeze(isCo.fValues(1,iNode,:));
                cond = obj.computeConformalMapCondition(isC);
                nodesC = connecC(:,iNode);
                nodesD = connecD(:,iNode);
                sC = sC + sparse(nodesD,nodesC,cond,nnodeD*nElemD,nnodesC);
            end
            obj.interpolator = sC;
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