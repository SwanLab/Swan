classdef EIFEMtesting < handle

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

        displacementFun
        LHS
        RHS
        scale
    end

    methods (Access = public)

        function obj = EIFEMtesting()
            close all
            obj.init();
            obj.createReferenceMesh();
           
            s.nsubdomains = obj.nSubdomains; %nx ny
            s.meshReference = obj.meshReference;         
            m = MeshCreatorFromRVE(s);
            [obj.meshDomain,meshSubDomain,interfaceConnec] = m.create();

            obj.createDisplacementFun();
            rawBC                  = obj.createRawBoundaryConditions();
            obj.boundaryConditions = obj.createBoundaryConditions(obj.meshDomain,rawBC);
            obj.quad               = Quadrature.set(obj.meshDomain.type);
            obj.quad.computeQuadrature('QUADRATIC');
            obj.createDomainMaterial();
            obj.computeStiffnessMatrix();
            obj.computeForces();
    
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
            obj.scale    = 'MACRO';
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

        function BC = createRawBoundaryConditions(obj)
            bM = obj.meshDomain.createBoundaryMesh();
            
            dirichletBc.boundaryId=1;
            dirichletBc.dof=[1,2];
            dirichletBc.value=[0,0];
            newmanBc.boundaryId=2;
            newmanBc.dof=[2];
            newmanBc.value=[10];

            [dirichlet,pointload] = obj.createBc(bM,dirichletBc,newmanBc);
            BC.dirichlet=dirichlet;
            BC.pointload=pointload;
%             obj.boundaryConditions = BC;
        end        

        function bc = createBoundaryConditions(obj,mesh,bcV)
            dim = obj.getFunDims();
            bcV.ndimf = dim.ndimf;
            bcV.ndofs = dim.ndofs;
            s.mesh  = mesh;
            s.scale = 'MACRO';
            s.bc    = {bcV};
            s.ndofs = dim.ndofs;
            bc = BoundaryConditions(s);
            bc.compute();
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

        function createDisplacementFun(obj)
            obj.displacementFun = P1Function.create(obj.meshDomain, obj.meshDomain.ndim);
        end

        function dim = getFunDims(obj)
            d.ndimf  = obj.displacementFun.ndimf;
            d.nnodes = size(obj.displacementFun.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.meshDomain.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function computeStiffnessMatrix(obj)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.meshDomain;
            s.fun      = obj.displacementFun;
            s.material = obj.material;
            lhs = LHSintegrator.create(s);
            obj.LHS = lhs.compute();
        end

        function computeForces(obj)
            s.type = 'Elastic';
            s.scale    = obj.scale;
            s.dim      = obj.getFunDims();
            s.BC       = obj.boundaryConditions;
            s.mesh     = obj.meshDomain;
            s.material = obj.material;
%             s.globalConnec = obj.displacementField.connec;
            s.globalConnec = obj.meshDomain.connec;
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
            R = RHSint.computeReactions(obj.LHS);
%             obj.variables.fext = rhs + R;
            obj.RHS = rhs;
        end

    end
end
