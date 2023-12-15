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
            dofsElements = obj.computeDofsElements(ndofsEdges);
            
            obj.dofs = [dofsVertices,dofsEdges,dofsElements];
            obj.ndofs = max(max(obj.dofs));
        end
        
        
        function computeCoordPriv(obj)
            nelem = size(obj.dofs,1);
            coor = zeros(obj.ndofs,2);
            
            if obj.order~=1
                for ielem = 1:nelem
                    nodesMesh = obj.mesh.coord(obj.mesh.connec(ielem,1:3),:);
                    coor(obj.dofs(ielem,:),:) = obj.computeNodesElement(nodesMesh);
                end
                obj.coord = coor;
            else
                obj.coord = obj.mesh.coord;
            end
        end
       
        
        function dofsVertices = computeDofsVertices(obj)
            m = obj.mesh;
            dofsVertices = m.connec;
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

                for iElem = 1:obj.mesh.nelem
                    for iNode = 1:obj.mesh.nnodeElem
                        ind = (iNode-1)*ndofEdge+1:iNode*ndofEdge;
                        if locPointEdge(iElem,iNode)~=iNode
                            ind = flip(ind);
                        end
                        dofsEdges(iElem,ind) = edges(iElem,iNode)*ndofEdge-(ndofEdge-1):edges(iElem,iNode)*ndofEdge;
                    end
                end
                dofsEdges = dofsEdges + m.nnodes;
            end
        end
        
        
        function dofsElements = computeDofsElements(obj,ndofsEdges)
            m = obj.mesh;
            polOrder = obj.order;
            if polOrder == 1
                dofsElements = [];
            else 
                ndofsElements = obj.computeNdofsElements(polOrder);
                dofsElements = zeros(m.nelem,ndofsElements);
                for i = 1:m.nelem
                    dofsElements(i,:) = (i-1)*ndofsElements+1:i*ndofsElements;
                end
                dofsElements = dofsElements + ndofsEdges;
            end
        end

        
        function ndofsElements = computeNdofsElements(obj,polOrder)
            ord = polOrder - 3;
            ndofsElements = 0;
            if ord == 0
                ndofsElements = ndofsElements + 1;
            elseif ord > 0
                ndofsElements = ndofsElements + obj.mesh.nnodeElem + obj.mesh.edges.nEdgeByElem*(ord-1);
                ndofsElements = ndofsElements + obj.computeNdofsElements(ord-3);
            end
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
        
        
        function coor = computeNodesElement(obj,coords)
            nodesBase = obj.interp.pos_nodes';          
            base = nodesBase(:,1:obj.mesh.nnodeElem);
            
            transfMatrix = (coords-coords(1,:))'/base;
            transfVect = coords(1,:)';
            coor = transfMatrix*nodesBase+transfVect;
            coor = coor';
        end
        
        
        function ord = convertOrder(~,order)
            switch order
                case 'P1'
                    ord = 1;
                case 'P2'
                    ord = 2;
                case 'P3'
                    ord = 3;
            end
        end

    end

end
