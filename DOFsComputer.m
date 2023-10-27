classdef DOFsComputer < handle
    
    properties (Access = public)
    end
    
    properties (Access = protected)
        dofs
        ndofs
        ndimf
        mesh
        order
    end
    
    methods (Access = public)
        
        function obj = DOFsComputer(cParams)
            obj.init(cParams);
        end
        
        function ndofs = computeNumberDofs(obj)
            dofsVertices = obj.computeDofsVertices();
            dofsEdges = obj.computeDofsEdges();
            dofsElements = obj.computeDofsElements(dofsEdges);
            
            obj.dofs = [dofsVertices,dofsEdges,dofsElements];
            obj.ndofs = max(max(obj.dofs));
            ndofs = obj.ndofs;
        end
        
    end
        
    methods (Access = protected)

        function init(obj,cParams)
            obj.ndimf = cParams.ndimf;
            obj.mesh = cParams.mesh;
            obj.order = obj.convertOrder(cParams.order);
        end
        
        function ord = convertOrder(~,order)
            switch order
                case 'LINEAR'
                    ord = 1;
                case 'QUADRATIC'
                    ord = 2;
                case 'CUBIC'
                    ord = 3;
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
                dofsEdges = dofsEdges + m.nnodes;
            end
        end
        
        function dofsElements = computeDofsElements(obj,dofsEdges)
            m = obj.mesh;
            polOrder = obj.order;
            if polOrder == 1
                dofsElements = [];
            else
                ord = polOrder-2;
                ndofsElements = 0;
                for i = 0:ord
                    ndofsElements = ndofsElements + i;
                end
                dofsElements = zeros(m.nelem,ndofsElements);
                for i = 1:m.nelem
                    dofsElements(i,:) = (i-1)*ndofsElements+1:i*ndofsElements;
                end
                dofsElements = dofsElements + max(max(dofsEdges));
            end
        end

    end

    methods (Abstract)
        
    end   

end
