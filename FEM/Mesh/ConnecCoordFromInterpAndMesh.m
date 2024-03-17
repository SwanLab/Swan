classdef ConnecCoordFromInterpAndMesh < handle

    properties (Access = public)
        coord
        connec
    end

    properties (Access = private)
        mesh
        interp
        order
    end

    methods (Access = public)

        function obj = ConnecCoordFromInterpAndMesh(cParams)
            obj.init(cParams)
        end

        function compute(obj)
            obj.computeDofs();
            obj.computeCoords();
%             obj.refineDofs();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh   = cParams.mesh;
            obj.interp = cParams.interpolation;
            obj.order   = cParams.order;
        end
        
        function computeDofs(obj)
            dofsVertices = obj.computeDofsVertices();
            dofsEdges = obj.computeDofsEdges();
            dofsElements = obj.computeDofsElements(dofsEdges);
            obj.connec = [dofsVertices,dofsEdges,dofsElements];
        end
        
        function refineDofs(obj)
            dofsDim =  obj.connec;
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
        
        function dofsVertices = computeDofsVertices(obj)
            dofsVertices = obj.mesh.connec;
        end
        
        function dofsEdges = computeDofsEdges(obj)
            obj.mesh.computeEdges();
            if obj.order == 1
                dofsEdges = [];
            else
                edges = obj.mesh.edges.edgesInElem;
                ndofEd = (obj.order-1);
                ndofsEdgeElem = ndofEd*obj.mesh.nnodeElem;
                dofsEdges = zeros(obj.mesh.nelem,ndofsEdgeElem);
                locPointEdge = squeeze(obj.mesh.edges.localNodeByEdgeByElem(:,:,1));

                for i = 1:obj.mesh.nelem
                    for j = 1:obj.mesh.nnodeElem
                        ind = (j-1)*ndofEd+1:j*ndofEd;
                        if (locPointEdge(i,j)~=j)
                            ind = flip(ind);
                        end
                        dofsEdges(i,ind) = edges(i,j)*ndofEd-(ndofEd-1):edges(i,j)*ndofEd;
                    end
                end
                dofsEdges = dofsEdges + obj.mesh.nnodes;
            end
        end
        
        function dofsElements = computeDofsElements(obj,dofsEdges)
            if obj.order == 1
                dofsElements = [];
            else
                ord = obj.order-2;
                ndofsElements = 0;
                for i = 0:ord
                    ndofsElements = ndofsElements + i;
                end
                dofsElements = zeros(obj.mesh.nelem,ndofsElements);
                for i = 1:obj.mesh.nelem
                    dofsElements(i,:) = (i-1)*ndofsElements+1:i*ndofsElements;
                end
                dofsElements = dofsElements + max(max(dofsEdges));
            end
        end
        
        function computeCoords(obj)
            nelem = size(obj.connec,1);
            ndofs = max(max(obj.connec));
            x = zeros(ndofs,1);
            y = zeros(ndofs,1);
            
            if obj.order~=1
                for ielem = 1:nelem
                    coor = obj.computeNodesElement(obj.mesh.coord(obj.mesh.connec(ielem,1:3),:));
                    x(obj.connec(ielem,:)) = coor(:,1);
                    y(obj.connec(ielem,:)) = coor(:,2);
                end
            
                obj.coord = [x,y];
            else
                obj.coord = obj.mesh.coord;
            end
        end
        
        function coor = computeNodesElement(obj,coords)
            base = obj.interp.pos_nodes;          
            c = base(1:3,:);
            M = (coords-coords(1,:))'/c';
            N = coords(1,:)';
            coor = (M*base'+N)';
        end

        function shapes = computeShapesInVariableNodes(obj,mesh)
            interpMesh = Interpolation.create(mesh.type,'LINEAR');
            nNodeMesh = interpMesh.nnode;
            nNodeVar  = obj.interp.nnode;
            shapes    = zeros(nNodeVar,nNodeMesh);
            nodesVar  = obj.interp.pos_nodes;
            for inodeVar = 1:obj.interp.nnode
                nodesPoints = nodesVar(inodeVar,:);
                shapes(inodeVar,:) = interpMesh.computeShapeFunctions(nodesPoints');
            end
        end

    end

    methods (Access = private, Static)

        function ind = findPointInList(node,xPoints)
            match = true(size(xPoints,1),1);
            for idime = 1:size(node,2)
%                 match = match & xPoints(:,idime) == node(idime);
                match = match & abs(xPoints(:,idime) - node(idime))<10e-8; % NOOOOOOO
            end
            ind = find(match);
        end

    end

end
