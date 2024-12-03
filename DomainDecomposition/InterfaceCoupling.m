classdef InterfaceCoupling < handle
    
    properties (Access = public)
        interfaceConnec
        interfaceConnecReshaped
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
            obj.computeCouplingConnec();
            obj.reshapeConecPerInterface();
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
            for jDom = 1:nY
                for iDom = 1:nX
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
            [GlNodeBd,ind]  = unique(GlNodeBd);
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
            obj.interfaceConnec=unique(sameNodeOrdered,'rows');
         end

         function reshapeConecPerInterface(obj)
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
            obj.interfaceConnecReshaped = intConec;
        end

       
    end
    
end