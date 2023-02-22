classdef SymmetricContMapCondition < handle
    
    properties (Access = private)
        isCoherent
        isCoherentFun
        orientationDisc
    end
    
    properties (Access = private)
       mesh
       orientationVector
    end
    
    methods (Access = public)
        
        function obj = SymmetricContMapCondition(cParams)
            obj.init(cParams);
        end
        
        function c = computeCondition(obj)
            obj.createOrientationDiscontinous();
            obj.isOrientationCoherent();
            c = obj.computeSymmetricCondition();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh              = cParams.mesh;
            obj.orientationVector = cParams.orientationVector;
        end
        
        function createOrientationDiscontinous(obj)
            f = obj.orientationVector;
            fD = f.project('P1D');
            obj.orientationDisc = fD.fValues;
        end

        function isOrientationCoherent(obj)
            s.mesh        = obj.mesh;
            s.orientation = obj.orientationDisc;
            c = CoherentOrientationSelector(s);
            isC = c.isOrientationCoherent();
            obj.isCoherent = isC;
            a.fValues = permute(double(isC), [3 2 1]);
            a.mesh = obj.mesh;
            obj.isCoherentFun = P1DiscontinuousFunction(a);
        end

        function sC = computeSymmetricCondition(obj)
            nnodeD    = obj.mesh.nnodeElem;
            nElemD    = obj.mesh.nelem;
            nnodesC   = obj.mesh.nnodes;
            connecC   = obj.mesh.connec;
            connecD = obj.computeDiscontinuousConnectivities();
            sC = sparse(nnodeD*nElemD,nnodesC);
            for iNode = 1:nnodeD
                isC  = obj.isCoherent(:,iNode);
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