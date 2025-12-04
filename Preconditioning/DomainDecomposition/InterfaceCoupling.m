classdef InterfaceCoupling < handle
    
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
        
        function obj = InterfaceCoupling(cParams)
            obj.init(cParams)
            
        end
        
        function compute(obj)
            obj.coordNodeBoundary();
            obj.computeCouplingConnec2();
         %   obj.reshapeConecPerInterface();
        end

      function intConec = reshapeConecPerInterface(obj)
          if (size(obj.interfaceConnec,2) == 2)
            globalConec = obj.interfaceConnec;
            nRnodes     = obj.meshReference.nnodes;
            nnode       = size(globalConec,1);
            nint        = obj.ninterfaces;
            gbInt        = 1;
            inode = 1;
            while inode < nnode
                nodeId = globalConec(inode,1);
                dom    = ceil(nodeId/nRnodes);
                row = ceil(dom/obj.nSubdomains(1));
                col = dom-(row-1)*obj.nSubdomains(1);
                bdMesh = obj.interfaceMeshSubDomain{row,col};
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
          else
              disp('Continuous mesh generated. However, case not suitable for domain decomposition yet')
              intConec = obj.interfaceConnec;
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
            nnodes        = obj.meshReference.nnodes;
            interfaceMesh = obj.interfaceMeshSubDomain();
            ndim          = interfaceMesh{1,1}{1,1}.mesh.ndim;
            ninterface    = obj.ninterfaces;
            coordBdGl     = zeros(1,ndim);
            GlNodeBd      = zeros(1,1);
                for jDom = 1:nX
                    for iDom = 1:nY
                        for iline=1:ninterface
                            bdcood    = interfaceMesh{jDom,iDom}{iline,1}.mesh.coord;
                            coordBdGl = [coordBdGl;bdcood];
                            %although it says global is in subdomain
                            %conecctivity
                            conecInter = interfaceMesh{jDom,iDom}{iline,1}.globalConnec;
                            nodeIntSub = reshape(unique(conecInter),[],1);
                            nodeIntGl  = nodeIntSub + nnodes*(nX*(jDom-1)+iDom-1);
                            GlNodeBd   = [GlNodeBd; nodeIntGl];
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
            ndim         = obj.meshReference.ndim;
            nBdNode      = length(obj.GlNodeBd);
            globalNode   = obj.GlNodeBd;
            coordBdGlAux = obj.coordBdGl;
            imaster=1;
            if ndim == 2
                coordAux = [coordBdGlAux zeros(nBdNode,1)];
            else
                coordAux = coordBdGlAux;
            end
            for iBdNode = 1:nBdNode
                NodeCoord = coordAux(iBdNode,:);
                tol = obj.tolSameNode;


                isSameNode = vecnorm(coordAux-NodeCoord,'Inf',2) - tol <= 0;

%                aux       = (coordAux(:,1)==NodeCoord(1) & coordAux(:,2)==NodeCoord(2) & coordAux(:,3)==NodeCoord(3));
%                 [~,ind]   = ismember(NodeCoord,coordAux,'Rows');
                %ind       = find(aux == 1);
                if sum(isSameNode)>1                   
                        sameNode = globalNode(isSameNode);
                        nsame   = length(sameNode);
                        sameNodeOrdered(imaster,1:nsame) = sort(sameNode);
                        imaster=imaster+1;                    
                end
            end
            obj.interfaceConnec= unique(sameNodeOrdered,'rows','stable');
%              obj.interfaceConnec= sameNodeOrdered;
         end
function computeCouplingConnec2(obj)
    ndim       = obj.meshReference.ndim;
    nBdNode    = numel(obj.GlNodeBd);
    globalNode = obj.GlNodeBd;
    coordBdGl  = obj.coordBdGl;
    tol        = obj.tolSameNode;

    % --- Normalize dimension (pad to 3D for uniformity)
    if ndim == 2
        coord = [coordBdGl zeros(nBdNode,1)];
    else
        coord = coordBdGl;
    end

    % --- Quantize coordinates by tolerance
    %     (nodes within tol map to same bin)
    quantized = round(coord / tol) * tol;

    % --- Find unique "bins" of coincident nodes
    [~, ~, ic] = unique(quantized, 'rows', 'stable');

    % --- Group nodes by unique coordinate bin
    % Preallocate max group size (usually small)
    nGroups = max(ic);
    maxGroupSize = accumarray(ic, 1, [nGroups,1], @max);
    maxSize = max(maxGroupSize);

    % Initialize interface connection table
    sameNodeOrdered = zeros(nGroups, maxSize);

    % Fill each group row with sorted node IDs
    for g = 1:nGroups
        sameNode = sort(globalNode(ic == g));
        sameNodeOrdered(g, 1:numel(sameNode)) = sameNode;
    end

    % --- Keep only rows with duplicates (same node sets)
    hasDuplicates = sum(sameNodeOrdered > 0, 2) > 1;
    obj.interfaceConnec = unique(sameNodeOrdered(hasDuplicates,:), 'rows', 'stable');
end

   

       
    end
    
end