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
        end
        
        
        function ndofs = getNumberDofs(obj)
            ndofs = obj.ndofs;
        end
        
        
        function dofs = getDofs(obj)
            dofs = obj.dofs;
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
            obj.ndofs = max(max(dofsFacelets));
        end

        function dofsFacelets = computeDofsFacelets(obj)
            if obj.mesh.ndim == 2
                dofsFacelets = obj.computeDofsEdges();
            elseif obj.mesh.ndim == 3
                dofsFacelets = obj.computeDofsFaces();
            end
        end
        
        function result = vectorUnion(~, v, n, fl)
            result = v + (0:n-1);
            result(fl,:) = flip(result(fl,:),2);
        end

        
        function dofsEdges = computeDofsEdges(obj)
            m = obj.mesh;
            m.computeEdges();
            edges = m.edges.edgesInElem;
            ndofEdge = 1;
            ndofsEdgeElem = obj.mesh.edges.nEdgeByElem;
            
            dofsEdges = zeros(obj.mesh.nelem,ndofsEdgeElem*ndofEdge);
            locPointEdge = squeeze(obj.mesh.edges.localNodeByEdgeByElem(:,:,1));
            locPointEdgeRef = obj.computeLocPointEdgeRef();

            for iEdge = 1:obj.mesh.edges.nEdgeByElem
                ind = (iEdge-1)*ndofEdge+1:iEdge*ndofEdge;
                fl = locPointEdge(:,iEdge)~=locPointEdgeRef(iEdge);
                v = edges(:,iEdge)*ndofEdge-(ndofEdge-1);
                dofsEdges(:,ind) = obj.vectorUnion(v, length(ind), fl);
            end
        end


        function dofsFaces = computeDofsFaces(obj)
            m = obj.mesh;
            m.computeFaces();
            dofsFaces = m.faces.facesInElem;
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
