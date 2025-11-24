classdef InterfaceCoupling3D < handle

    properties (Access = public)
        interfaceConnec
    end

    properties (Access = private)
        meshReference
        nSubdomains
        interfaceMeshSubDomain
        ninterfaces
        tolSameNode
    end

    properties (Access = private)
        coordBdGl
        GlNodeBd
    end

    methods (Access = public)

        function obj = InterfaceCoupling3D(cParams)
            obj.init(cParams)

        end

        function compute(obj)
            obj.coordNodeBoundary();
            obj.computeCouplingConnec();
            %   obj.reshapeConecPerInterface();
        end

        function intConec = reshapeConecPerInterface(obj)
            globalConec = obj.interfaceConnec;
            nRnodes     = obj.meshReference.nnodes;
            nnode       = size(globalConec,1);
            nint        = obj.ninterfaces;
            gbInt        = 1;
            inode = 1;
            while inode < nnode
                nodeId = globalConec(inode,1);
                dom    = ceil(nodeId/nRnodes);
                kInd   = ceil(dom/(obj.nSubdomains(1)*obj.nSubdomains(2)));
                pageD  = dom-(kInd-1)*obj.nSubdomains(1)*obj.nSubdomains(2);
                row = ceil(pageD/obj.nSubdomains(1));
                col = pageD-(row-1)*obj.nSubdomains(1);
                bdMesh = obj.interfaceMeshSubDomain{row,col,kInd};
                for iInt = 1:nint
                    isNode = bdMesh{iInt}.globalConnec == nodeId - nRnodes*(dom-1);
                    if sum(sum(isNode)) > 0
                        nnodebd = bdMesh{iInt}.mesh.nnodes;
                        intConec{gbInt} = globalConec(inode:inode+nnodebd-1,:);
                        inode = inode + nnodebd;
                        gbInt = gbInt+1;
                        break
                    end
                end
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.meshReference = cParams.meshReference;
            obj.nSubdomains   = cParams.nSubdomains;
            obj.interfaceMeshSubDomain = cParams.interfaceMeshSubDomain;
            obj.ninterfaces   = cParams.ninterfaces;
            obj.tolSameNode   = cParams.tolSameNode;
        end

        function coordNodeBoundary(obj)
            nX            = obj.nSubdomains(1);
            nY            = obj.nSubdomains(2);
            nZ            = obj.nSubdomains(3);
            nnodes        = obj.meshReference.nnodes;
            interfaceMesh = obj.interfaceMeshSubDomain();
            ndim          = interfaceMesh{1,1}{1,1}.mesh.ndim;
            ninterface    = obj.ninterfaces;
            coordBdGl     = zeros(1,ndim);
            GlNodeBd      = zeros(1,1);
            for kDom = 1:nZ
                for jDom = 1:nY
                    for iDom = 1:nX
                        for iline=1:ninterface
                            bdcood    = interfaceMesh{jDom,iDom,kDom}{iline,1}.mesh.coord;
                            coordBdGl = [coordBdGl;bdcood];
                            %although it says global is in subdomain
                            %conecctivity
                            conecInter = interfaceMesh{jDom,iDom,kDom}{iline,1}.globalConnec;
                            nodeIntSub = reshape(unique(conecInter),[],1);
                            nodeIntGl  = nodeIntSub + nnodes*(nX*nY*(kDom-1)+nX*(jDom-1)+iDom-1);
                            GlNodeBd   = [GlNodeBd; nodeIntGl];
                        end
                    end
                end
            end
            [GlNodeBd,ind]  = unique(GlNodeBd,'stable');
            coordBdGl       = coordBdGl(ind,:);
            obj.coordBdGl   = coordBdGl(2:end,:);
            %             subDomNode = subDomNode(2:end,:);
            obj.GlNodeBd  = GlNodeBd(2:end,:);
        end

        function computeCouplingConnec(obj)
            % I leave this chunk of code commented to recall wha i was
            % doing. Uncommented is chatgpt acceleration.!!!
            
%             ndim         = obj.meshReference.ndim;
%             nBdNode      = length(obj.GlNodeBd);
%             globalNode   = obj.GlNodeBd;
%             coordBdGlAux = obj.coordBdGl;
%             imaster=1;
%             if ndim == 2
%                 coordAux = [coordBdGlAux zeros(nBdNode,1)];
%             else
%                 coordAux = coordBdGlAux;
%             end
%             for iBdNode = 1:nBdNode
%                 NodeCoord = coordAux(iBdNode,:);
%                 tol = obj.tolSameNode;
% 
% 
%                 isSameNode = vecnorm(coordAux-NodeCoord,'Inf',2) - tol <= 0;
% 
%                 %                aux       = (coordAux(:,1)==NodeCoord(1) & coordAux(:,2)==NodeCoord(2) & coordAux(:,3)==NodeCoord(3));
%                 %                 [~,ind]   = ismember(NodeCoord,coordAux,'Rows');
%                 %ind       = find(aux == 1);
%                 if sum(isSameNode)>1
%                     sameNode = globalNode(isSameNode);
%                     nsame   = length(sameNode);
%                     sameNodeOrdered(imaster,1:nsame) = sort(sameNode);
%                     imaster=imaster+1;
%                 end
%             end
%             obj.interfaceConnec= unique(sameNodeOrdered,'rows','stable');

            ndim         = obj.meshReference.ndim;
            nBdNode      = length(obj.GlNodeBd);
            globalNode   = obj.GlNodeBd;
            coordBdGlAux = obj.coordBdGl;
            tol          = obj.tolSameNode;

            % Ensure 3D coordinates for uniformity
            if ndim == 2
                coordAux = [coordBdGlAux zeros(nBdNode,1)];
            else
                coordAux = coordBdGlAux;
            end

            % Use rangesearch with k-d tree (fast and memory efficient)
            Mdl = createns(coordAux, 'Distance', 'chebychev');  % or 'euclidean'
            idxList = rangesearch(Mdl, coordAux, tol);

            % Create edge list for nodes within tolerance (excluding self)
            edges = [];
            for i = 1:nBdNode
                neighbors = idxList{i};
                neighbors(neighbors == i) = [];  % Remove self
                if ~isempty(neighbors)
                    edges = [edges; repmat(i, length(neighbors), 1), neighbors(:)];
                end
            end

            % Group connected nodes using graph components
            if isempty(edges)
                obj.interfaceConnec = [];  % No duplicates
                return;
            end

            G = graph(edges(:,1), edges(:,2));
            bins = conncomp(G);  % Connected component labels

            % Build ordered groupings
            nGroups = max(bins);
            maxGroupSize = max(histcounts(bins,1:nGroups+1));
            sameNodeOrdered = zeros(nGroups, maxGroupSize);

            for i = 1:nGroups
                nodes = globalNode(bins == i);
                sameNodeOrdered(i, 1:length(nodes)) = sort(nodes);
            end

            % Final output
            rowSums = sum(sameNodeOrdered > 0, 2);
            sameNodeOrderedFiltered = sameNodeOrdered(rowSums > 1, :);
            obj.interfaceConnec = unique(sameNodeOrderedFiltered, 'rows', 'stable');
            %              obj.interfaceConnec= sameNodeOrdered;

        end




    end

end