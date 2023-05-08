classdef DomainDecompositionManager < handle

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (Access = private)
        meshDomain
        boundaryConditions
        material
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
            obj.createDomainMesh();
            obj.createBoundaryConditions();
            obj.createDomainMaterial();
            obj.solveDomainProblem();
        end
        
        function solveDomainProblem(obj)
            s.mesh     = obj.meshDomain;
            s.bc       = obj.boundaryConditions;
            s.material = obj.material;
            s.type = 'ELASTIC';
            s.scale = 'MACRO';
            s.dim = '2D';
            fem = FEM.create(s);
            fem.solve();
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

       function createDomainMesh(obj)
            %aquestes dos primeres en una classe q es digui SubdomainNodeRelator o algo aixi i q et torni el "MasterSlave"...
            % pero q no es digui masterSlave, millor nom (llenguatge %
            % inclusiu) ok?  no, master i slave no esta ben vist, millor
            % driver i reference o no se pero un altre nom
            %el tema de la esclvitud!

            [coordBdGl,subNodeDof,GlNodeBd] = obj.coordNodeBoundary();
            obj.createMasterSlaveConnec(coordBdGl,subNodeDof,GlNodeBd);

            % Una altre clase per crear Domain mesh, ok?
            connecGlob = obj.createGlobalConnec();
            coordGlob  = obj.createGlobalCoord();
            connecGlob = obj.updateGlobalConnec(connecGlob);
            coordGlob  = obj.updateGlobalCoord(coordGlob);
            s.connec   = connecGlob;
            s.coord    = coordGlob;
            m    = Mesh(s);
            m.plot()
            obj.meshDomain = m;
            figure
            m.plot
        end

        function createDomainMaterial(obj)
            ngaus = 1;        
            m = obj.meshDomain;                      
            obj.material = obj.createMaterial(m,ngaus);
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
                        GlNodeBd   = [GlNodeBd; nodeIntGl];
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

        function createMasterSlaveConnec(obj,coordBdGl,GlNodeBd)
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

        function createBoundaryConditions(obj)
            bM = obj.meshDomain.createBoundaryMesh();
            
            dirichletBc.boundaryId=1;
            dirichletBc.dof=[1,2];
            dirichletBc.value=[0,0];
            newmanBc.boundaryId=2;
            newmanBc.dof=[2];
            newmanBc.value=[10];

            [dirichlet,pointload] = obj.createBc(bM,dirichletBc,newmanBc);
            bc.dirichlet=dirichlet;
            bc.pointload=pointload;   
            obj.boundaryConditions = bc;
        end        


        function [dirichlet,pointload] = createBc(obj,boundaryMesh,dirchletBc,newmanBc)
            dirichlet = obj.createBondaryCondition(boundaryMesh,dirchletBc);
            pointload = obj.createBondaryCondition(boundaryMesh,newmanBc);
        end

        function cond = createBondaryCondition(obj,bM,condition)
            nbound = length(condition.boundaryId);
            cond = zeros(1,3);
            for ibound=1:nbound
                ncond  = length(condition.dof(nbound,:));
                nodeId= unique(bM{condition.boundaryId(ibound)}.globalConnec);
                nbd   = length(nodeId);
                for icond=1:ncond
                    bdcond= [nodeId, repmat(condition.dof(icond),[nbd,1]), repmat(condition.value(icond),[nbd,1])];
                    cond=[cond;bdcond];
                end
            end
            cond = cond(2:end,:);
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
