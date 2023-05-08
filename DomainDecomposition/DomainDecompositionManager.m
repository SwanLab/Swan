classdef DomainDecompositionManager < handle

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (Access = private)
        meshReference
        interfaceMeshReference
        meshSubDomain
        nSubdomains
        ninterfaces
        interfaceMeshSubDomain
        globalMeshConnecSubDomain
        interfaceMeshConnecSubDomain
        subDomainContact
        MasterSlaveConnec
        cornerNodes
    end

    methods (Access = public)

        function obj = DomainDecompositionManager()
            close all
            obj.init();
            obj.createReferenceMesh();
            obj.obtainCornerNodes();
            obj.createSubDomainMeshes();
            obj.createInterfaceSubDomainMeshes();
            %             obj.findContact();
            [coordBdGl,subNodeDof,GlNodeBd] = obj.coordNodeBoundary();
            obj.createMasterSlaveConnec(coordBdGl,subNodeDof,GlNodeBd);
            connecGlob = obj.createGlobalConnec();
            coordGlob  = obj.createGlobalCoord();
            connecGlob = obj.updateGlobalConnec(connecGlob);
            coordGlob  = obj.updateGlobalCoord(coordGlob);
            s.connec=connecGlob;
            s.coord = coordGlob;
            mFem=Mesh(s);
            mFemBd=mFem.createBoundaryMesh();
            figure
            mFem.plot

            dirichletBc.boundaryId=1;
            dirichletBc.dof=[1,2];
            dirichletBc.value=[0,0];

            newmanBc.boundaryId=2;
            newmanBc.dof=[2];
            newmanBc.value=[10];

            [dirichlet,pointload] = obj.createBc(mFemBd,dirichletBc,newmanBc);

            cParams.mesh=mFem;
            cParams.bc.dirichlet=dirichlet;
            cParams.bc.pointload=pointload;
            ngaus=1;
            cParams.material = createMaterial(obj,mFem,ngaus);
            cParams.scale = 'MACRO';
            cParams.dim = '2D';
            femp = FEM(cParams);


            %             obj.mesh        = cParams.mesh;
            %             obj.material    = cParams.material;
            %             obj.scale       = cParams.scale;
            %             obj.pdim        = cParams.dim;
            %             obj.inputBC     = cParams.bc;
            %             obj.computeTotalConnec();
            %             connecGlob = concatenate(obj);
            %             obj.computeGlobalConnec();
            % obj.createSubDomainStiffnessMatrices();
            % obj.computeGlobalStiffnessMatrix();
        end

    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains = [3 3]; %nx ny
        end

        function createReferenceMesh(obj)
            filename='lattice_ex1';
            a.fileName=filename;
            femD = FemDataContainer(a);
            mS = femD.mesh;
            bS = mS.createBoundaryMesh();
            obj.meshReference = mS;
            obj.interfaceMeshReference = bS;
            obj.ninterfaces=length(bS);
        end

        function L = computeReferenceMeshLength(obj)
            coord = obj.meshReference.coord;
            Lx = max(coord(:,1));
            Ly = max(coord(:,2));
            L = [Lx Ly];
        end

        function createSubDomainMeshes(obj)
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            figure(2)
            for jDom = 1:nY
                for iDom = 1:nX
                    coordIJ = obj.computeSubdomainCoords(jDom,iDom);
                    mIJ     = obj.createSubdomainMesh(coordIJ);
                    mIJ.plot();
                    mD{jDom,iDom} = mIJ;
                    hold on
                end
            end
            obj.meshSubDomain = mD;
        end

        function m = createSubdomainMesh(obj,coord)
            connec0  = obj.meshReference.connec;
            s.coord  = coord;
            s.connec = connec0;
            m = Mesh(s);
        end


        function coord = computeSubdomainCoords(obj,jDom,iDom)
            coord0 = obj.meshReference.coord;
            L  = obj.computeReferenceMeshLength();
            Lx = L(1);
            Ly = L(2);
            coord(:,1) = coord0(:,1)+Lx*(iDom-1);
            coord(:,2) = coord0(:,2)+Ly*(jDom-1);
        end

        function createInterfaceSubDomainMeshes(obj)
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            figure
            for jDom = 1:nY
                for iDom = 1:nX
                    bIJ = obj.meshSubDomain{jDom,iDom}.createBoundaryMesh();
                    bD{jDom,iDom} = bIJ;
                    hold on
                    for iline=1:length(bIJ)
                        bIJ{iline}.mesh.plot();
                    end
                end
            end
            obj.interfaceMeshSubDomain = bD;
        end


        function  computeTotalConnec(obj)
            %             connec0 = obj.meshReference.connec;
            %             nnodes  = obj.meshReference.nnodes;
            nX            = obj.nSubdomains(1);
            nY            = obj.nSubdomains(2);
            %             connecGlob    = obj.createGlobalConnec();
            contact       = obj.subDomainContact();
            interfaceMesh = obj.interfaceMeshSubDomain();
            for jDom = 1:nY
                for iDom = 1:nX
                    if ~isempty(contact.index)

                    end
                    %                     dConnec = connec0 + nnodes*(nX*(jDom-1)+iDom-1);
                    %                     gcoord =
                    gD{jDom,iDom} = dConnec;
                end
            end
            obj.globalMeshConnec = gD;
        end

        function findContact(obj)
            %only those bottom and left
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            %             ninterf=obj.ninterfaces;
            for jDom = 1:nY
                for iDom = 1:nX
                    con=1;
                    if jDom-1>0
                        contact(jDom,iDom).index(con,:)= [jDom-1 iDom];
                        contact(jDom,iDom).line(con) = 3;
                        con=con+1;
                    end
                    if iDom-1>0
                        contact(jDom,iDom).index(con,:)= [jDom iDom-1];
                        contact(jDom,iDom).line(con) = 1;
                    end
                end
            end
            obj.subDomainContact = contact;
        end

        function connecGlob = createGlobalConnec(obj)
            % we simply create a connectivity matrix as if no nodes are
            % shared, just summing nnodes. We need that to create the
            % global connectivity matrix that considers shared nodes.
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            nelem=obj.meshReference.nelem;
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
        end

        function  coordGlob=createGlobalCoord(obj)
            % we simply create the coordinate matrix as if no nodes are
            % shared, just summing nnodes.
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            %             nelem=obj.meshReference.nelem;
            %             nnodeElem=obj.meshReference.nnodeElem;
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
        end

        function [coordBdGl,subDomNode,GlNodeBd] = coordNodeBoundary(obj)
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            %             nelem=obj.meshReference.nelem;
            %             nnodeElem=obj.meshReference.nnodeElem;
            nnodes  = obj.meshReference.nnodes;
            %             connec0 =obj.meshReference.connec;
            %             dConnec = connec0 + nnodes*(nX*(jDom-1)+iDom-1);
            %             connecGlob=zeros(nX*nY*nelem,nnodeElem);
            interfaceMesh = obj.interfaceMeshSubDomain();
            ndim          = interfaceMesh{1,1}{1,1}.mesh.ndim;
            ninterface    = obj.ninterfaces();
            coordBdGl     = zeros(1,ndim);
            subDomNode     = zeros(1,1);
            GlNodeBd         = zeros(1,1);
            for jDom = 1:nY
                for iDom = 1:nX
                    for iline=1:ninterface
                        bdcood    = interfaceMesh{jDom,iDom}{iline,1}.mesh.coord;
                        coordBdGl = [coordBdGl;bdcood];
                        %although it says global is in subdomain
                        %conecctivity
                        conecInter  = interfaceMesh{jDom,iDom}{iline,1}.globalConnec;
                        nodeIntSub    = unique(conecInter);
                        nodeIntGl  = nodeIntSub + nnodes*(nX*(jDom-1)+iDom-1);
                        subDomNode = [subDomNode;nodeIntSub];
                        GlNodeBd     = [GlNodeBd; nodeIntGl];
                        %                     indLinear= nX*(jDom-1)+iDom;
                        %                     rowIn=(indLinear-1)*nelem+1;
                        %                     rowEnd=indLinear*nelem;
                        %                     connecGlob(rowIn:rowEnd,:)=connec0+nnodes*(indLinear-1);
                    end
                end
            end
            coordBdGl = coordBdGl(2:end,:);
            subDomNode = subDomNode(2:end,:);
            GlNodeBd   = GlNodeBd(2:end,:);
        end

        function createMasterSlaveConnec(obj,coordBdGl,subNodeDof,GlNodeBd)
            ndim      = obj.meshReference.ndim;
            nBdNode   = length(GlNodeBd);
            GlNodeAux = GlNodeBd;
            imaster=1;
            if ndim == 2
                coordAux = [coordBdGl zeros(nBdNode,1)];
            else
                coordAux = coordBdGl;
            end
            for iBdNode = 1:nBdNode
                NodeCoord = coordAux(iBdNode,:);
                aux       = (coordAux(:,1)==NodeCoord(1) & coordAux(:,2)==NodeCoord(2) & coordAux(:,3)==NodeCoord(3));
                ind       = find(aux == 1);
                if length(ind)>1
                    sameNode_aux = GlNodeAux(ind);
                    sameNode(imaster,:) = sort(sameNode_aux);
                    imaster=imaster+1;
                end
            end
            obj.MasterSlaveConnec=unique(sameNode,'rows');
        end

        function  connecGlob = updateGlobalConnec(obj,connecGlob)
            connecAux=connecGlob;
            mSconnec=obj.MasterSlaveConnec;
            nmaster=size(mSconnec,1);
            nslave = size(mSconnec,2);
            for imaster=1:nmaster
                for islave=2:nslave
                    connecGlob(connecGlob==mSconnec(imaster,islave))=mSconnec(imaster,1);
                    connecGlob(connecGlob>mSconnec(imaster,islave))=connecGlob(connecGlob>mSconnec(imaster,islave))-1;
                    mSconnec(mSconnec>mSconnec(imaster,islave))=mSconnec(mSconnec>mSconnec(imaster,islave))-1;
                end
            end
        end

        function  coordGlob = updateGlobalCoord(obj,coordGlob)
            coordGlob  = unique(coordGlob,'rows','stable');
        end

        function obtainCornerNodes(obj)
            ninterface = obj.ninterfaces;
            corner = zeros(ninterface,1);
            icorner=1;
            for i=1:ninterface-1
                nodesi=obj.interfaceMeshReference{i,1}.globalConnec;
                for j=i+1:ninterface
                    nodesj=obj.interfaceMeshReference{j,1}.globalConnec;
                    nodes=intersect(nodesi,nodesj);
                    if ~isempty(nodes)
                        corner(icorner)= nodes;
                        icorner=icorner+1;
                    end
                end
            end
            obj.cornerNodes = corner;
        end

        %         function obtainFaceNodes(obj)
        %             ninterface = obj.ninterfaces;
        %             corner = obj.cornerNodes;
        %             icorner=1;
        %             if corner(1)>0
        %               for i=1:ninterface
        %                 nodesi=obj.interfaceMeshReference{i,1}.globalConnec;
        % %                 nodesi(nodesi==)
        %               end
        %             else
        %
        %             end
        %
        %
        %             dist2face = abs(x - xLim);
        %             isInFace = dist2face < 1e-13;
        %             nodesInFace = obj.allNodes(isInFace);
        %         end

        function [dirichlet,pointload] = createBc(obj,boundaryMesh,dirchletBc,newmanBc)
            dirichlet = obj.createDirichletBc(boundaryMesh,dirchletBc);
            pointload = obj.createNewmanBc(boundaryMesh,newmanBc);
        end

        function dirichlet = createDirichletBc(obj,boundaryMesh,dirichletBc)
            nbound = length(dirichletBc.boundaryId);
            dirichlet = zeros(1,3);
            for ibound=1:nbound
                ncond  = length(dirichletBc.dof(nbound,:));
                nodeId= unique(boundaryMesh{dirichletBc.boundaryId(ibound)}.globalConnec);
                nbd   = length(nodeId);
                for icond=1:ncond
                    bdcond= [nodeId, repmat(dirichletBc.dof(icond),[nbd,1]), repmat(dirichletBc.value(icond),[nbd,1])];
                    dirichlet=[dirichlet;bdcond];
                end
            end
            dirichlet = dirichlet(2:end,:);
        end

        function pointload = createNewmanBc(obj,boundaryMesh,newmanBc)
            nbound = length(newmanBc.boundaryId);
            pointload = zeros(1,3);
            for ibound=1:nbound
                ncond  = length(newmanBc.dof(nbound,:));
                nodeId= unique(boundaryMesh{newmanBc.boundaryId(ibound)}.globalConnec);
                nbd   = length(nodeId);
                for icond=1:ncond
                    bdcond= [nodeId, repmat(newmanBc.dof(icond),[nbd,1]), repmat(newmanBc.value(icond),[nbd,1])];
                    pointload=[pointload;bdcond];
                end
            end
            pointload = pointload(2:end,:);
        end

        function material = createMaterial(obj,mesh,ngaus)
            I = ones(mesh.nelem,ngaus);
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = mesh.nelem;
            s.mesh  = mesh;
            s.kappa = .9107*I;
            s.mu    = .3446*I;
            mat = Material.create(s);
            mat.compute(s);
            material = mat;
        end
    end
end
