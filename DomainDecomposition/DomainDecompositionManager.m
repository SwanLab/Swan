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
        cornerNodes
        quad
    end

    methods (Access = public)

        function obj = DomainDecompositionManager()
            close all
            obj.init();
            obj.createReferenceMesh();
           
            s.nsubdomains = obj.nSubdomains; %nx ny
            s.meshReference = obj.meshReference;         
            m = MeshCreatorFromRVE(s);
            [meshDomain,meshSubDomain,interfaceConnec] = m.create();
    
            cMesh = createCoarseMesh(obj);
            s.nsubdomains = obj.nSubdomains; %nx ny
            s.meshReference = cMesh;         
            mRVECoarse = MeshCreatorFromRVE(s);
            [meshDomainCoarse,meshSubDomainCoarse,interfaceConnecCoarse] = mRVECoarse.create();

            filename = 'testEIFEM1.mat';
            RVE = TrainedRVE(filename);
            s1.RVE = RVE;
            s1.mesh = meshDomainCoarse;
            eifem = EIFEM(s1);
            
% %             obj.obtainCornerNodes();
% %             fineMesh = MeshFromRVE
%             obj.createSubDomainMeshes();
%             obj.createInterfaceSubDomainMeshes();
%             obj.createDomainMesh();
            obj.createBoundaryConditions();
            obj.quad = Quadrature.set(obj.meshDomain.type);
            obj.quad.computeQuadrature('QUADRATIC');
            obj.createDomainMaterial();
            
%             s.referenceMesh = obj.referenceMesh;
%             mC = MeshCreatorFromSubmeshes();
%             obj.meshDomain = mC.mesh;

%             preconditioner = obj.createPreconditioner(mC.submeshes);

            obj.solveDomainProblem();



        end
        
        function solveDomainProblem(obj)
            s.mesh     = obj.meshDomain;
            s.bc       = obj.boundaryConditions;
            s.material = obj.material;
            s.type     = 'ELASTIC';
            s.scale    = 'MACRO';
            s.dim      = '2D';
            s.solverTyp = 'ITERATIVE';
            s.iterativeSolverTyp = 'PCG';
            s.preconditionerType = 'EIFEM';
            s.tol = 1e-6;
            
            fem        = FEM.create(s);
            fem.solve();
        end

    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains = [2 2]; %nx ny
        end

        function createReferenceMesh(obj)
%             filename   = 'lattice_ex1';
%             a.fileName = filename;
%             femD       = FemDataContainer(a);
%             mS         = femD.mesh;
%             bS         = mS.createBoundaryMesh();
             % Generate coordinates
            x1 = linspace(0,1,5);
            x2 = linspace(0,1,5);
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectivities
            [F,coord] = mesh2tri(xv,yv,zeros(size(xv)),'x');

            s.coord    = coord(:,1:2);
            s.connec   = F;
            mS         = Mesh(s);
            bS         = mS.createBoundaryMesh();

            obj.meshReference = mS;
            obj.interfaceMeshReference = bS;
            obj.ninterfaces = length(bS);
        end

        function cMesh = createCoarseMesh(obj)
            xmax = max(obj.meshReference.coord(:,1));
            xmin = min(obj.meshReference.coord(:,1));
            ymax = max(obj.meshReference.coord(:,2));
            ymin = min(obj.meshReference.coord(:,2));
            coord(1,1) = xmax;
            coord(1,2) = ymin;
            coord(2,1) = xmax;
            coord(2,2) = ymax;
            coord(3,1) = xmin;
            coord(3,2) = ymax;
            coord(4,1) = xmin;
            coord(4,2) = ymin;
            connec = [1 2 3 4];
            s.coord = coord;
            s.connec = connec;
            cMesh = Mesh(s);
        end

       function createDomainMesh(obj)
            s.nSubdomains   = obj.nSubdomains;
            s.meshReference = obj.meshReference;
            s.interfaceMeshSubDomain = obj.interfaceMeshSubDomain();
            s.ninterfaces   = obj.ninterfaces;
            s.meshSubDomain = obj.meshSubDomain;

            coupling = InterfaceCoupling(s);
            coupling.compute();
            s.interfaceConnec = coupling.interfaceConnec;

            DMesh = DomainMeshComputer(s);
            DMesh.compute();
            obj.meshDomain = DMesh.domainMesh;
        end

        function createDomainMaterial(obj)
%             ngaus = 1;        
            m = obj.meshDomain;                      
            obj.material = obj.createMaterial(m);
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

        function obtainCornerNodes(obj)
            ninterface = obj.ninterfaces;
%             corner = zeros(ninterface,1);
            icorner=1;
%             for i=1:ninterface-1
%                 nodesi=obj.interfaceMeshReference{i,1}.globalConnec;
%                 for j=i+1:ninterface
%                     nodesj=obj.interfaceMeshReference{j,1}.globalConnec;
%                     nodes=intersect(nodesi,nodesj);
%                     if ~isempty(nodes)
%                         corner(icorner)= nodes;
%                         icorner=icorner+1;
%                     end
%                 end
%             end
            for i=1:ninterface
                aux = CellNodesDescriptor(obj.interfaceMeshReference{i}.mesh.coord);
                corner = unique(aux.cornerNodes);
                

                [coordmax]  = max(obj.interfaceMeshReference{i}.mesh.coord);
                [coordmin]  = max(obj.interfaceMeshReference{i}.mesh.coord);

                corner(icorner) = indmax;
                icorner         = icorner+1;
                
                corner(icorner) = indmin;
                icorner         = icorner+1;
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

        function material = createMaterial(obj,mesh)
            I = ones(mesh.nelem,obj.quad.ngaus);
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
