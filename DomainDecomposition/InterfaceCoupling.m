classdef InterfaceCoupling < handle
    
    properties (Access = public)
        interfaceConnec
    end
    
    properties (Access = private)
        meshReference
        nSubdomains
        interfaceMeshSubDomain
        ninterfaces
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
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.meshReference = cParams.meshReference;
            obj.nSubdomains   = cParams.nSubdomains;
            obj.interfaceMeshSubDomain = cParams.interfaceMeshSubDomain;
            obj.ninterfaces   = cParams.ninterfaces;
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
            GlNodeAux    = obj.GlNodeBd;
            coordBdGlAux = obj.coordBdGl;
            imaster=1;
            if ndim == 2
                coordAux = [coordBdGlAux zeros(nBdNode,1)];
            else
                coordAux = coordBdGlAux;
            end
            for iBdNode = 1:nBdNode
                NodeCoord = coordAux(iBdNode,:);
%                 aux       = ((abs(coordAux(:,1)-NodeCoord(1))<=1e-14) & (abs(coordAux(:,2)-NodeCoord(2))<=1e-14) & (abs(coordAux(:,3)-NodeCoord(3))<=1e-14));
                aux       = (coordAux(:,1)==NodeCoord(1) & coordAux(:,2)==NodeCoord(2) & coordAux(:,3)==NodeCoord(3));
%                 [~,ind]   = ismember(NodeCoord,coordAux,'Rows');
                ind       = find(aux == 1);
                if length(ind)>1                   
                        sameNode_aux = GlNodeAux(ind);
                        nsame = length(sameNode_aux);
                        sameNode(imaster,1:nsame) = sort(sameNode_aux);
                        imaster=imaster+1;                    
                end
            end
            obj.interfaceConnec=unique(sameNode,'rows');
         end
        
    end
    
end