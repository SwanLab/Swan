classdef DomainMeshComputer < handle

    properties (Access = public)
        domainMesh
        localGlobalConnec
    end

    properties (Access = private)
        meshReference
        nSubdomains
        interfaceMeshSubDomain
        ninterfaces
        meshSubDomain
        interfaceConnec
        
    end

    properties (Access = private)
        connecGlob
        coordGlob
        updtConnecGlob
        updtCoordGlob
    end

    methods (Access = public)

        function obj = DomainMeshComputer(cParams)
            obj.init(cParams)

        end

        function compute(obj)
            obj.createGlobalConnec();
            obj.createGlobalCoord();
            obj.updateGlobalConnec();
            obj.updateGlobalCoord();
            obj.computeLocalGlobalConnec()
            s.coord        = obj.updtCoordGlob;
            s.connec       = obj.updtConnecGlob;
            obj.domainMesh = Mesh.create(s); 
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.meshReference   = cParams.meshReference;
            obj.nSubdomains     = cParams.nSubdomains;
            obj.interfaceMeshSubDomain = cParams.interfaceMeshSubDomain;
            obj.ninterfaces     = cParams.ninterfaces;
            obj.meshSubDomain   = cParams.meshSubDomain;
            obj.interfaceConnec = cParams.interfaceConnec;
        end

        function createGlobalConnec(obj)
            % we simply create a connectivity matrix as if no nodes are
            % shared, just summing nnodes. We need that to create the
            % global connectivity matrix that considers shared nodes.
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            nelem = obj.meshReference.nelem;
            nnodeElem=obj.meshReference.nnodeElem;
            nnodes  = obj.meshReference.nnodes;
            connec0 =obj.meshReference.connec;
            %             dConnec = connec0 + nnodes*(nX*(jDom-1)+iDom-1);
            connecGlob=zeros(nX*nY*nelem,nnodeElem);
            for jDom = 1:nY
                for iDom = 1:nX
                    indLinear= nX*(jDom-1)+iDom;
                    rowIn=(indLinear-1)*nelem+1;
                    rowEnd=indLinear*nelem;
                    connecGlob(rowIn:rowEnd,:)=connec0+nnodes*(indLinear-1);
                end
            end
            obj.connecGlob=connecGlob;
        end

        function  createGlobalCoord(obj)
            % we simply create the coordinate matrix as if no nodes are
            % shared, just summing nnodes.
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            ndim     = obj.meshReference.ndim;
            nnodes  = obj.meshReference.nnodes;
            meshsd = obj.meshSubDomain;
            %             dConnec = connec0 + nnodes*(nX*(jDom-1)+iDom-1);
            coordGlob=zeros(nX*nY*nnodes,ndim);
            for jDom = 1:nY
                for iDom = 1:nX
                    indLinear= nX*(jDom-1)+iDom;
                    rowIn=(indLinear-1)*nnodes+1;
                    rowEnd=indLinear*nnodes;
                    coordGlob(rowIn:rowEnd,:)=meshsd{jDom,iDom}.coord;
                end
            end
            obj.coordGlob=coordGlob;
        end

        function   updateGlobalConnec(obj)
            updtConnecGlob = obj.connecGlob;
            rCconnec       = obj.interfaceConnec;
            nref           = size(rCconnec,1);
%             ncopies        = size(rCconnec,2);
            for iref=1:nref
                aux = rCconnec(iref,:);
                aux = aux(aux>0);
                for icopy=2:length(aux)
                    updtConnecGlob(updtConnecGlob==aux(icopy))  = aux(1);
                    updtConnecGlob(updtConnecGlob>aux(icopy))   = updtConnecGlob(updtConnecGlob>aux(icopy))-1;
                    aux(aux>aux(icopy)) =  aux(aux>aux(icopy))-1;
                    rCconnec(rCconnec>aux(icopy)) = rCconnec(rCconnec>aux(icopy))-1;
                end
            end
            obj.updtConnecGlob=updtConnecGlob;
        end

        function  updateGlobalCoord(obj)
            obj.updtCoordGlob  = unique(obj.coordGlob,'rows','stable');
        end

        function  computeLocalGlobalConnec(obj)
            nodeG  = reshape(obj.updtConnecGlob',[],1);
            nodeL  = reshape(obj.connecGlob',[],1);
            nnode = obj.meshReference.nnodes;
            localGlobalConnec(:,1) = nodeG;
            localGlobalConnec(:,2) = nodeL;
            localGlobalConnec      = unique(localGlobalConnec,'rows','stable');    
            localGlobalConnec      = permute(reshape(localGlobalConnec',2,nnode,[]),[2,1,3]);
           
            for dom = 1:size(localGlobalConnec,3)
                localGlobalConnec(:,2,dom) = localGlobalConnec(:,2,dom)-(dom-1)*nnode;
            end
             obj.localGlobalConnec  = localGlobalConnec;
        end
    end
end