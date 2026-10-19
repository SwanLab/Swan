classdef DomainMeshComputer < handle

    properties (Access = public)
        domainMesh
        domainMeshDisc
        localGlobalConnec
    end

    properties (Access = private)
        meshReference
        nSubdomains
        interfaceMeshSubDomain
        ninterfaces
        meshSubDomain
        interfaceConnec
        tolSameNode
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
            obj.connecGlob = obj.createGlobalConnec(obj.meshReference);
            obj.createGlobalCoord();
            obj.updateGlobalConnec();
            obj.updateGlobalCoord();
            obj.computeLocalGlobalConnec()
            s.coord        = obj.updtCoordGlob;
            s.connec       = obj.updtConnecGlob;
            obj.domainMesh = Mesh.create(s);
            ss.connec = obj.connecGlob;
            ss.coord  = obj.coordGlob;
            obj.domainMeshDisc = Mesh.create(ss);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.tolSameNode            = cParams.tolSameNode;
            obj.meshReference          = cParams.meshReference;
            obj.nSubdomains            = cParams.nSubdomains;
            obj.interfaceMeshSubDomain = cParams.interfaceMeshSubDomain;
            obj.ninterfaces            = cParams.ninterfaces;
            obj.meshSubDomain          = cParams.meshSubDomain;
            obj.interfaceConnec        = cParams.interfaceConnec;
        end

        function gConnec = createGlobalConnec(obj,m)
            % we simply create a connectivity matrix as if no nodes are
            % shared, just summing nnodes. We need that to create the
            % global connectivity matrix that considers shared nodes.
            nX         = obj.nSubdomains(1);
            nY         = obj.nSubdomains(2);
            nelem      = m.nelem;
            nnodeElem  = m.nnodeElem;
            nnodes     = m.nnodes;
            connec0    = m.connec;
            gConnec    = zeros(nX*nY*nelem,nnodeElem);
            for jDom = 1:nY
                for iDom = 1:nX
                    indLinear               = nX*(jDom-1)+iDom;
                    rowIn                   = (indLinear-1)*nelem+1;
                    rowEnd                  = indLinear*nelem;
                    gConnec(rowIn:rowEnd,:) = connec0+nnodes*(indLinear-1);
                end
            end
        end

        function  createGlobalCoord(obj)
            % we simply create the coordinate matrix as if no nodes are
            % shared, just summing nnodes.
            nX        = obj.nSubdomains(1);
            nY        = obj.nSubdomains(2);
            ndim      = obj.meshReference.ndim;
            nnodes    = obj.meshReference.nnodes;
            meshsd    = obj.meshSubDomain;
            coordGlob = zeros(nX*nY*nnodes,ndim);
            for jDom = 1:nY
                for iDom = 1:nX
                    indLinear= nX*(jDom-1)+iDom;
                    rowIn=(indLinear-1)*nnodes+1;
                    rowEnd=indLinear*nnodes;
                    coordGlob(rowIn:rowEnd,:)=meshsd{jDom,iDom}.coord;
                end
            end
            obj.coordGlob = coordGlob;
        end
        
        function updateGlobalConnec(obj)
            updtConnecGlob = obj.connecGlob;
            rCconnec       = obj.interfaceConnec;
            maxID          = max(max(updtConnecGlob(:)), max(rCconnec(:)));
            map            = 1:maxID;
            nref           = size(rCconnec, 1);
            for iref = 1:nref
                aux = rCconnec(iref, :);
                aux = aux(aux > 0);
                if length(aux) > 1
                    master_node = aux(1);
                    slave_nodes = aux(2:end);
                    map(map(slave_nodes)) = map(master_node);
                end
            end
            updtConnecGlob      = map(updtConnecGlob);
            unique_map_targets  = unique(map);
            renum_map           = zeros(1, maxID);
            renum_map(unique_map_targets) = 1:length(unique_map_targets);

            updtConnecGlob     = renum_map(updtConnecGlob);
            obj.updtConnecGlob = updtConnecGlob;
        end

        function  updateGlobalCoord(obj)
            tol   = obj.tolSameNode;
            nodes = obj.connecGlob(:);
            nodes = unique(nodes);
            for i = 1:size(obj.interfaceConnec,2)
                [~,ind]           = ismember(obj.interfaceConnec(:,i),nodes);
                nodes(ind(ind>0)) = obj.interfaceConnec(ind>0,1);
            end
            nodes             = unique(nodes,'stable');
            uniqueVals        = obj.coordGlob(nodes,:);
            obj.updtCoordGlob = uniqueVals;
        end


        function  computeLocalGlobalConnec(obj)
            nodeG  = reshape(obj.updtConnecGlob',[],1);
            nodeL  = reshape(obj.connecGlob',[],1);
            nnode  = obj.meshReference.nnodes;
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