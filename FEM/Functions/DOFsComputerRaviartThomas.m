classdef DOFsComputerRaviartThomas < handle
    
    properties (Access = public)
    end
    
    properties (Access = private)
        dofs
        ndofs
        ndimf
        mesh
        interp
        dofsOrientation
    end
    
    methods (Access = public)
        
        function obj = DOFsComputerRaviartThomas(cParams)
            obj.init(cParams);
        end
        
        
        function computeDofs(obj)
            obj.computeDofsPriv();
            obj.refineDofs();
            obj.computeDofsOrientation();
        end
        
        
        function ndofs = getNumberDofs(obj)
            ndofs = obj.ndofs;
        end
        
        
        function dofs = getDofs(obj)
            dofs = obj.dofs;
        end

        function dofsOri = getDofsOrientation(obj)
            dofsOri = obj.dofsOrientation;
        end
        
    end
        
    methods (Access = private)

        function init(obj,cParams)
            obj.ndimf = cParams.ndimf;
            obj.mesh = cParams.mesh;
            
            if any(contains(fieldnames(cParams),'interpolation'))
                obj.interp = cParams.interpolation;
            end
        end
        
        
        function computeDofsPriv(obj)
            dofsFacelets = obj.computeDofsFacelets();
            obj.dofs = dofsFacelets;
        end

        function dofsFacelets = computeDofsFacelets(obj)
            if obj.mesh.ndim == 2
                dofsFacelets = obj.computeDofsEdges();
            elseif obj.mesh.ndim == 3
                dofsFacelets = obj.computeDofsFaces();
            end
        end
        
        function dofsEdges = computeDofsEdges(obj)
            m = obj.mesh;
            m.computeEdges();
            edges = m.edges.edgesInElem;
            ndofEdge = 1;
            ndofsEdgeElem = obj.mesh.edges.nEdgeByElem;
            
            dofsEdges = zeros(obj.mesh.nelem,ndofsEdgeElem);
            locPointEdge = squeeze(obj.mesh.edges.localNodeByEdgeByElem(:,:,1));
            locPointEdgeRef = obj.computeLocPointEdgeRef();

            for iElem = 1:obj.mesh.nelem
                for iEdge = 1:obj.mesh.edges.nEdgeByElem
                    ind = (iEdge-1)*ndofEdge+1:iEdge*ndofEdge;
                    if locPointEdge(iElem,iEdge)~=locPointEdgeRef(iEdge)
                        ind = flip(ind);
                    end
                    dofsEdges(iElem,ind) = edges(iElem,iEdge)*ndofEdge-(ndofEdge-1):edges(iElem,iEdge)*ndofEdge;
                end
            end
        end


        function dofsFaces = computeDofsFaces(obj,ndofsEdges)
            m = obj.mesh;
            m.computeFaces();
            dofsFaces = m.faces.facesInElem;
        end


        function computeDofsOrientation(obj)

            m = obj.mesh;
            m.computeEdges();
            ndofsEdgeElem = obj.mesh.edges.nEdgeByElem;

            dofsOri = ones(obj.mesh.nelem,ndofsEdgeElem);
            locPointEdge = squeeze(obj.mesh.edges.localNodeByEdgeByElem(:,:,1));
            locPointEdgeRef = obj.computeLocPointEdgeRef();

            for iElem = 1:obj.mesh.nelem
                for iEdge = 1:obj.mesh.edges.nEdgeByElem
                    if locPointEdge(iElem,iEdge)~=locPointEdgeRef(iEdge)
                        dofsOri(iElem,iEdge) = -1;
                    end
                end
            end

            obj.dofsOrientation = dofsOri;

        end


        function refineDofs(obj)
            dofsDim =  obj.dofs;
            nDimf  = obj.ndimf;
            nNode  = size(dofsDim, 2);
            nDofsE = nNode*nDimf;
            dofsElem  = zeros(nDofsE,size(dofsDim,1));
            for iNode = 1:nNode
                for iUnkn = 1:nDimf
                    idofElem   = nDimf*(iNode - 1) + iUnkn;
                    globalNode = dofsDim(:,iNode);
                    idofGlobal = nDimf*(globalNode - 1) + iUnkn;
                    dofsElem(idofElem,:) = idofGlobal;
                end
            end
            obj.dofs = dofsElem';
            obj.ndofs = max(max(obj.dofs));
        end
        
        function loc = computeLocPointEdgeRef(obj)
            type = obj.mesh.type;
            switch type
                case 'LINE'
                    loc = [1 2];
                case 'TRIANGLE'
                    loc = [1 2 3];
                case 'QUAD'
                    loc = [1 2 3 4];
                case 'TETRAHEDRA'
                    loc = [1 1 1 2 2 3];
                case 'HEXAHEDRA'
                    loc = [1 4 1 2 2 3 3 4 5 8 6 7];
            end
        end
        

    end

end
