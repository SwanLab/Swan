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
        continuityType
    end
    
    methods (Access = public)
        
        function obj = DOFsComputer(cParams)
            obj.init(cParams);
        end
        
        
        function computeDofs(obj)
            obj.computeDofsPriv();
            obj.obtainDofsForVectorField();
        end
        
        
        function computeCoord(obj)
            if isempty(obj.dofs)
                obj.computeDofs();
            end

            if  isa(obj.mesh,'Mesh') 
                obj.computeCoordPriv(obj.dofs);
            end
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

            if isequal(cParams.order,'P1D')
                obj.continuityType = 'Disc';
            else
                obj.continuityType = 'Cont';
            end
            
            if any(contains(fieldnames(cParams),'interpolation'))
                obj.interp = cParams.interpolation;
            end
        end
        
        
        function computeDofsPriv(obj)
            dofsVertices = obj.computeDofsVertices();
            
            dofsEdges = obj.computeDofsEdges();
            ndofsEdges = max(max(dofsEdges));
            
            if ~strcmp(obj.mesh.type,'LINE')
                    dofsFaces = obj.computeDofsFaces(ndofsEdges);
                ndofsFaces = max(max(dofsFaces));
            
                dofsElements = obj.computeDofsElements(ndofsFaces);
            else
                dofsFaces = [];
                dofsElements = [];
            end

            obj.dofs = [dofsVertices,dofsEdges,dofsFaces,dofsElements];
        end
        
        
        function computeCoordPriv(obj,dofs)
            switch obj.continuityType
                case 'Cont'
                    ndofsE = size(dofs,2);
                    if obj.order~=1
                        coor   = zeros(obj.ndofs,obj.mesh.ndim);
                        sAF      = obj.computefHandlePosition();
                        sAF.mesh = obj.mesh;
                        sAF.ndimf = obj.ndimf;
                        func     = AnalyticalFunction(sAF);
                        c = func.evaluate(obj.interp.pos_nodes');
                        for idof = 1:ndofsE
                            for idim = 1:obj.mesh.ndim
                                nodes    = dofs(:,idof);
                                iNode = ceil(idof/obj.ndimf);
                                coor(nodes,idim) =  c(idim,iNode,:);
                            end
                        end
                        obj.coord = coor;
                    else
                        coor   = zeros(obj.ndofs,obj.mesh.ndim);
                        for idim = 1:obj.mesh.ndim
                            nodes    = unique(obj.mesh.connec);
                            coorDofs = repmat(obj.mesh.coord(nodes,idim)',obj.ndimf,1);
                            dofs     = obj.computeNodesToDofs(nodes);
                            coor(dofs,idim) = coorDofs(:);
                        end
                        obj.coord = coor;
                    end
                case 'Disc'
                    coordC = obj.mesh.coord;
                    connec = obj.mesh.connec;
                    dofsC  = reshape(connec',1,[]);
                    coorD  = zeros(obj.ndofs,obj.mesh.ndim);
                    for idim = 1:obj.mesh.ndim
                        coorI = repmat(coordC(dofsC,idim),1,obj.ndimf)';
                        coorD(1:obj.ndofs,idim) = reshape(coorI,[],1);
                    end
                    obj.coord = coorD;
            end

            end

        function dofs = computeNodesToDofs(obj,nodes)
            for iDimf = 1:obj.ndimf
                dofs(iDimf,:) =(nodes-1)*obj.ndimf+iDimf;
            end
            dofs = dofs(:);
        end

        
        function dofsVertices = computeDofsVertices(obj)
            switch obj.continuityType
                case 'Cont'
                    if obj.order == 0
                        dofsVertices = (1:obj.mesh.nelem)';
                    else
                        m = obj.mesh;
                        dofsVertices = m.connec;
                    end
                case 'Disc'
                nDofs = obj.mesh.nnodeElem*obj.mesh.nelem;                    
                nodes = 1:nDofs;
                dofsVertices = reshape(nodes,obj.mesh.nnodeElem,obj.mesh.nelem)';        
            end
        end


        function result = vectorUnion(~, v, n, fl)

            result = v + (0:n-1);
            result(fl,:) = flip(result(fl,:),2);
        end

        
        function dofsEdges = computeDofsEdges(obj)
            if obj.order <= 1
                dofsEdges = [];
            else
                m = obj.mesh;
                m.computeEdges();
                edges = m.edges.edgesInElem;
                ndofEdge = obj.order-1;
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

                dofsEdges = dofsEdges + m.nnodes;
            end
        end
        
        
        function dofsFaces = computeDofsFaces(obj,ndofsEdges)
            m = obj.mesh;
            polOrder = obj.order;
            if polOrder <= 1
                dofsFaces = [];
            else 
                m.computeFaces();
                faces = m.faces.facesInElem;
                ndofFace = obj.computeNdofsFaces(polOrder);
                ndofsFaceElem = ndofFace*obj.mesh.faces.nFaceByElem;
                
                dofsFaces = zeros(obj.mesh.nelem,ndofsFaceElem);
                for iElem = 1:obj.mesh.nelem
                    for iFace = 1:obj.mesh.faces.nFaceByElem
                        ind = (iFace-1)*ndofFace+1:iFace*ndofFace;
                        dofsFaces(iElem,ind) = faces(iElem,iFace)*ndofFace-(ndofFace-1):faces(iElem,iFace)*ndofFace;
                    end
                end
                
                dofsFaces = dofsFaces + ndofsEdges;
            end
        end
        
        function ndofsFaces = computeNdofsFaces(obj,polOrder)
            switch obj.mesh.type
                case {'TRIANGLE','TETRAHEDRA'}
                    d = 3;
                case {'QUAD','HEXAHEDRA'}
                    d = 2;
            end
            
            ord = polOrder - d;
            ndofsFaces = 0;
            if ord == 0
                ndofsFaces = ndofsFaces + 1;
            elseif ord > 0
                ndofsFaces = ndofsFaces + obj.mesh.faces.nNodeByFace + obj.mesh.faces.nEdgeByFace*(ord-1);
                ndofsFaces = ndofsFaces + obj.computeNdofsFaces(ord-d);
            end
        end
        
        
        function dofsElements = computeDofsElements(obj,ndofsFaces)
            m = obj.mesh;
            polOrder = obj.order; 
            isNot3D = strcmp(m.type,'TRIANGLE') || strcmp(m.type,'LINE') || strcmp(m.type,'QUAD'); 
            if polOrder <= 1 || isNot3D
                dofsElements = [];
            else 
                ndofElement = obj.computeNdofsElements(polOrder);
                dofsElements = zeros(obj.mesh.nelem,ndofElement);
                for iElem = 1:obj.mesh.nelem
                    dofsElements(iElem,:) = (iElem-1)*ndofElement+1:iElem*ndofElement;
                end
                
                dofsElements = dofsElements + ndofsFaces;
            end
        end
        
        function ndofsElements = computeNdofsElements(obj,polOrder)
            switch obj.mesh.type
                case 'TETRAHEDRA'
                    d = 4;
                case 'HEXAHEDRA'
                    d = 2;
            end
            
            ord = polOrder - d;
            ndofsElements = 0;
            if ord == 0
                ndofsElements = ndofsElements + 1;
            elseif ord > 0
                r.coord = obj.mesh.coord(obj.mesh.connec(1,:),:);
                r.connec = 1:size(r.coord,1);
                s.mesh = Mesh.create(r);
                s.order   = ord;
                s.ndimf   = obj.ndimf;
                c = DOFsComputer(s);
                c.computeDofs();
                ndofsElements = ndofsElements + c.getNumberDofs();
            end
        end
        
        
        function obtainDofsForVectorField(obj)
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
        
        
        function ord = convertOrder(~,order)
            if ischar(order)
                switch order
                    case 'P0'
                        ord = 0;
                    case {'P1','P1D'}
                        ord = 1;
                    case 'P2'
                        ord = 2;
                    case 'P3'
                        ord = 3;
                end
            else
                ord = order;
            end
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
        
        function fh = computefHandlePosition(obj)
            ndim = obj.mesh.ndim;
            switch ndim
                case 1
                    fh.fHandle = @(x) x(1,:,:);
                    fh.ndimf = 1;
                case 2
                    fh.fHandle = @(x) [x(1,:,:);x(2,:,:)];
                    fh.ndimf = 2;
                case 3
                    fh.fHandle = @(x) [x(1,:,:);x(2,:,:);x(3,:,:)];
                    fh.ndimf = 3;
            end
        end

    end

end
