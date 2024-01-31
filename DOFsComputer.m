classdef DOFsComputer < handle
    
    properties (Access = public)
    end
    
    properties (Access = private)
        dofs
        ndofs
        ndimf
        mesh
        order
        coord
        interp
    end
    
    methods (Access = public)
        
        function obj = DOFsComputer(cParams)
            obj.init(cParams);
        end
        
        
        function computeDofs(obj)
            obj.computeDofsPriv();
            obj.refineDofs();
        end
        
        
        function computeCoord(obj)
            if isempty(obj.dofs)
                obj.computeDofs;
            end
            obj.computeCoordPriv();
        end
        
        
        function ndofs = getNumberDofs(obj)
            ndofs = obj.ndofs;
        end
        
        
        function dofs = getDofs(obj)
            dofs = obj.dofs;
        end
        
        
        function coord = getCoord(obj)
            coord = obj.coord;
        end
        
    end
        
    methods (Access = private)

        function init(obj,cParams)
            obj.ndimf = cParams.ndimf;
            obj.mesh = cParams.mesh;
            obj.order = obj.convertOrder(cParams.order);
            
            if any(contains(fieldnames(cParams),'interpolation'))
                obj.interp = cParams.interpolation;
            end
        end
        
        
        function computeDofsPriv(obj)
            dofsVertices = obj.computeDofsVertices();
            dofsEdges = obj.computeDofsEdges();
            ndofsEdges = max(max(dofsEdges));
            dofsFaces = obj.computeDofsFaces(ndofsEdges);
            
            obj.dofs = [dofsVertices,dofsEdges,dofsFaces];
            obj.ndofs = max(max(obj.dofs));
        end
        
        
        function computeCoordPriv(obj)
            nelem = size(obj.dofs,1);
            coor = zeros(obj.ndofs,obj.mesh.ndim);
            
            if obj.order~=1
                
                sAF.fHandle = obj.computefHandlePosition();
                sAF.ndimf   = 2;
                sAF.mesh    = obj.mesh;
                func = AnalyticalFunction(sAF);
                c = func.evaluate(obj.interp.pos_nodes');
                c = reshape(c,obj.interp.nnode,obj.mesh.ndim,obj.mesh.nelem);
                
                for ielem = 1:nelem
                    coor(obj.dofs(ielem,:),:) = c(:,:,ielem);
                end
                
                obj.coord = coor;
            else
                obj.coord = obj.mesh.coord;
            end
        end
       
        
        function dofsVertices = computeDofsVertices(obj)
            if obj.order == 0
                dofsVertices = (1:obj.mesh.nelem)';
            else
                m = obj.mesh;
                dofsVertices = m.connec;
            end
        end
        
        
        function dofsEdges = computeDofsEdges(obj)
            m = obj.mesh;
            m.computeEdges();
            
            if obj.order == 1
                dofsEdges = [];
            else
                edges = m.edges.edgesInElem;
                ndofEdge = obj.order-1;
                ndofsEdgeElem = ndofEdge*obj.mesh.nnodeElem;
                
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
                dofsEdges = dofsEdges + m.nnodes;
            end
        end
        
        
        function dofsElements = computeDofsFaces(obj,ndofsEdges)
            m = obj.mesh;
            polOrder = obj.order;
            if polOrder == 1
                dofsElements = [];
            else 
                ndofsElements = obj.computeNdofsFaces(polOrder);
                dofsElements = zeros(m.nelem,ndofsElements);
                for i = 1:m.nelem
                    dofsElements(i,:) = (i-1)*ndofsElements+1:i*ndofsElements;
                end
                dofsElements = dofsElements + ndofsEdges;
            end
        end

        
        function ndofsElements = computeNdofsFaces(obj,polOrder)
            switch obj.mesh.type
                case 'TRIANGLE'
                    d = 3;
                    nfaces = 1;
                case 'QUAD'
                    d = 2;
                    nfaces = 1;
                case 'TETRAHEDRA'
                    d = 3;
                    nfaces = 4;
            end
            
            ord = polOrder - d;
            ndofsElements = 0;
            if ord == 0
                ndofsElements = ndofsElements + 1;
            elseif ord > 0
                ndofsElements = ndofsElements + obj.mesh.nnodeElem + obj.mesh.edges.nEdgeByElem*(ord-1);
                ndofsElements = ndofsElements + obj.computeNdofsElements(ord-d);
            end
            
            ndofsElements = ndofsElements*nfaces;
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
        end
        
        
        function ord = convertOrder(~,order)
            switch order
                case 'P0'
                    ord = 0;
                case 'P1'
                    ord = 1;
                case 'P2'
                    ord = 2;
                case 'P3'
                    ord = 3;
            end
        end
        
        
        function loc = computeLocPointEdgeRef(obj)
            type = obj.mesh.type;
            switch type
                case 'TRIANGLE'
                    loc = [1 2 3];
                case 'QUAD'
                    loc = [1 2 3 4];
                case 'TETRAHEDRA'
                    loc = [1 1 1 2 2 3];
            end
        end
        
        function fh = computefHandlePosition(obj)
            ndim = obj.mesh.ndim;
            switch ndim
                case 1
                    fh = @(x) x(1,:,:);
                case 2
                    fh = @(x) [x(1,:,:),x(2,:,:)];
                case 3
                    fh = @(x) [x(1,:,:),x(2,:,:),x(3,:,:)];
            end
        end

    end

end
