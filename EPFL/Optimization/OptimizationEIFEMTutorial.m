classdef OptimizationEIFEMTutorial < handle

    properties (Access = private)
        mesh
        filter
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        volume
        cost
        constraint
        dualVariable
        optimizer

        referenceMesh 
        boundaryRefMesh
        coarseMesh
        subdomainMeshes
        nSubdomains
        tolSameNode
        r
        xmin
        xmax
        ymin
        ymax
        Nr
        Ntheta
        x0
        y0
        BC
        sDir
        bcApplier
        fileNameEIFEM
        EIFEM
        solverType
        interfaceConnec 
        localGlobal     
        iCR             
        discMesh        
        EIFEMprecontitioner
        volumeTarget
        primalUpdater
    end

    methods (Access = public)

        function obj = OptimizationEIFEMTutorial()
            obj.init()
            obj.createMesh();
            obj.createCoarseMesh();
            obj.createDesignVariable();
%             obj.createFilter();
%             obj.createMaterialInterpolator();
            [LHSr,RHSr]= obj.createElasticProblem();
            obj.createEIFEMPreconditioner(RHSr);
%             obj.createComplianceFromConstiutive();
%             obj.createCompliance();
            obj.createComplianceRadius();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createDualVariable();
            obj.createPrimalUpdater()
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
            obj.nSubdomains = [28,14]; %50 15
            rmin = 1e-6;
            obj.r = rmin*ones(obj.nSubdomains)'; 
%             obj.r= (1e-6 - 1e-6) * rand(obj.nSubdomains(2),obj.nSubdomains(1)) + 1e-6;
            obj.xmax=1; obj.xmin=-1; obj.ymax = 1; obj.ymin=-1; 
            obj.Nr = 7; obj.Ntheta = 14; % for circle/square
%             obj.Nr = 10; obj.Ntheta = 10;% for lattice
            obj.x0 = 0; obj.y0=0;
            obj.tolSameNode = 1e-10;
%             obj.fileNameEIFEM = './EPFL/parametrizedEIFEMLagrange20_der2_lattice.mat';
%             obj.fileNameEIFEM = './EPFL/parametrizedEIFEMLagrange20_der2_092.mat';
            obj.fileNameEIFEM = './EPFL/parametrizedEIFEMLagrange20_der2_092.mat';
            obj.solverType = 'REDUCED';
            obj.volumeTarget = 0.6; %0.7
        end

        function createMesh(obj)
            mSbd = obj.createSubDomainMeshes();
            bS = mSbd{1,1}.createBoundaryMesh();
            [mD,mSb,iC,lG,iCR,discMesh] = obj.createMeshDomainJoiner(mSbd);
            obj.mesh = mD;
            obj.subdomainMeshes = mSb;
            obj.interfaceConnec = iC;
            obj.localGlobal     = lG;
            obj.iCR             = iCR;
            obj.discMesh        = discMesh;
        end

         function  mSbd = createSubDomainMeshes(obj)
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            Lx = obj.xmax-obj.xmin;
            Ly = obj.ymax-obj.ymin;
            for jDom = 1:nY
                for iDom = 1:nX
                    refMesh = mesh_rectangle_via_triangles(obj.r(jDom,iDom),obj.xmax,obj.xmin,obj.ymax,obj.ymin,obj.Nr,obj.Ntheta,obj.x0,obj.y0);
%                     refMesh = mesh_square_X_solid(1,obj.r(jDom,iDom),obj.Nr,obj.Ntheta);
                    coord0 = refMesh.coord;
                    s.coord(:,1) = coord0(:,1)+Lx*(iDom-1);
                    s.coord(:,2) = coord0(:,2)+Ly*(jDom-1);
                    s.connec = refMesh.connec;
                    mIJ     = Mesh.create(s);
                    %                     plot(mIJ)
                    %                     hold on;
                    mSbd{jDom,iDom} = mIJ;
                end
            end
            obj.referenceMesh = mSbd{1,1};
            obj.boundaryRefMesh = obj.referenceMesh.createBoundaryMesh();
        end

         function [mD,mSb,iC,lG,iCR,discMesh] = createMeshDomainJoiner(obj,mSbd)
            s.nsubdomains   = obj.nSubdomains;
            s.meshReference = obj.referenceMesh;
            s.tolSameNode = obj.tolSameNode;
            s.meshSbd     = mSbd;
            %             m = MeshCreatorFromRVE.create(s);
            m = MeshJoiner(s);
            [mD,mSb,iC,~,lG,iCR,discMesh] = m.create();
         end

          function createCoarseMesh(obj)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = obj.createReferenceCoarseMesh();
            s.tolSameNode   = obj.tolSameNode;
            mRVECoarse      = MeshCreatorFromRVE.create(s);
            [obj.coarseMesh,~,~] = mRVECoarse.create();  
        end

        function cMesh = createReferenceCoarseMesh(obj)
            coord(1,1) = obj.xmin;  coord(1,2) = obj.ymin;
            coord(2,1) = obj.xmax;  coord(2,2) = obj.ymin;
            coord(3,1) = obj.xmax;  coord(3,2) = obj.ymax;
            coord(4,1) = obj.xmin;  coord(4,2) = obj.ymax;
            connec = [2 3 4 1];
%                                 connec = [1 2 3 4];
            s.coord = coord;
            s.connec = connec;
            cMesh = Mesh.create(s);
        end

        function createDesignVariable(obj)
            s.ndimf    = 1;
            s.fValues  = obj.r(:);
            s.mesh     = obj.coarseMesh;
            s.order    = 'P0';
            s.fun      = LagrangianFunction(s);
            s.type     = 'Radius';
            s.plotting = true;
            s.nSubdomains = obj.nSubdomains;
            s.Nr       = obj.Nr;
            s.Ntheta   = obj.Ntheta;
            s.xmax     = obj.xmax;
            s.xmin     = obj.xmin;
            s.ymax     = obj.ymax;
            s.ymin     = obj.ymin;
            s.x0       = obj.x0;
            s.y0       = obj.y0;
            s.discMesh = obj.discMesh;
            radius     = DesignVariable.create(s);
            obj.designVariable = radius;
        end

        function createFilter(obj)
            s.filterType = 'P1';
            s.mesh  = obj.coarseMesh;
            s.trial = LagrangianFunction.create(obj.coarseMesh,1,'P0');
            s.test = LagrangianFunction.create(obj.coarseMesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

        function createMaterialInterpolator(obj)
            E0 = 1e-3;
            nu0 = 1/3;
            ndim = obj.mesh.ndim;
            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);


            E1 = 1;
            nu1 = 1/3;
            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;
        end

%         function m = createMaterial(obj)
%             x = obj.designVariable.fun;           
%             s.type                 = 'DensityBased';
%             s.density              = x;
%             s.materialInterpolator = obj.materialInterpolator;
%             s.dim                  = '2D';
%             s.mesh                 = obj.mesh;
%             m = Material.create(s);
%         end

        function material = createMaterial(obj)
            E  = 1;
            nu = 1/3;
            %              Epstr  = E/(1-nu^2);
            %             nupstr = nu/(1-nu);
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = obj.mesh.ndim;
            s.young   = ConstantFunction.create(E,obj.mesh);
            s.poisson = ConstantFunction.create(nu,obj.mesh);
            tensor    = Material.create(s);
            material  = tensor;
        end

%         function createElasticProblem(obj)
%             s.mesh = obj.mesh;
%             s.scale = 'MACRO';
%             s.material = obj.createMaterial();
%             s.dim = '2D';
%             s.boundaryConditions = obj.createBoundaryConditions();
%             s.interpolationType = 'LINEAR';
%             s.solverType = 'REDUCED';
%             s.solverMode = 'DISP';
%             s.solverCase = 'DIRECT';
%             fem = ElasticProblem(s);
%             obj.physicalProblem = fem;
%         end

         function [LHSr,RHSr] = createElasticProblem(obj)
            obj.BC = obj.createBoundaryConditions(); 
            obj.createBCapplier();
            material = obj.createMaterial();
            u = LagrangianFunction.create(obj.mesh,obj.mesh.ndim,'P1');            
            [lhs,LHSr] = obj.computeStiffnessMatrix(obj.mesh,u,material);
            [RHS,RHSr]       = obj.computeForces(lhs,u);
         end

           function [LHS,LHSr] = computeStiffnessMatrix(obj,mesh,dispFun,C)
            LHS = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u))),dispFun,dispFun,mesh,'Domain',2);
            LHSr = obj.bcApplier.fullToReducedMatrixDirichlet(LHS);
           end

             function [RHS, RHSr] =  computeForces(obj,stiffness,u)
            bc  = obj.BC;
            t   = bc.tractionFun;
            rhs = zeros(u.nDofs,1);
            if ~isempty(t)
                for i = 1:numel(t)
                    rhsi = t(i).computeRHS(u);
                    rhs  = rhs + rhsi;
                end
            end
            if strcmp(obj.solverType,'REDUCED')
                bc      = obj.BC;
                dirich  = bc.dirichlet_dofs;
                dirichV = bc.dirichlet_vals;
                if ~isempty(dirich)
                    R = -stiffness(:,dirich)*dirichV;
                else
                    R = zeros(sum(obj.uFun.nDofs(:)),1);
                end
                rhs = rhs+R;
            end
            RHS  = rhs;
            RHSr = obj.bcApplier.fullToReducedVectorDirichlet(rhs);
        end

        function createBCapplier(obj)
            s.mesh                  = obj.mesh;
            s.boundaryConditions    = obj.BC;
            obj.bcApplier           = BCApplier(s);
        end

        function Meifem = createEIFEMPreconditioner(obj,RHSr)
            mR            = obj.referenceMesh;
            s.RVE         = TrainedRVE(obj.fileNameEIFEM);
            s.mesh        = obj.coarseMesh;
            s.DirCond     = obj.sDir;
            s.nSubdomains = obj.nSubdomains;
            s.mu          = obj.r;
            s.meshRef     = obj.discMesh;
%             s.Fext        = RHSr;
            eifem         = EIFEMnonPeriodic(s);
%             obj.EIFEMprecontitioner = eifem;

            iC  = obj.interfaceConnec;
            lG  = obj.localGlobal;
            iCR = obj.iCR; 
            bS  = obj.boundaryRefMesh;
            ss.ddDofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,mR,iCR);
            ss.EIFEMsolver  = eifem;
            ss.bcApplier    = obj.bcApplier;
            ss.dMesh        = obj.discMesh;
            ss.type         = 'EIFEM';
            ss.Fext         = RHSr;
            eP     = Preconditioner.create(ss);
            obj.EIFEMprecontitioner = eP;
%             obj.EIFEMprecontitioner = @(r) eP.apply(r);
        end

         function d = createDomainDecompositionDofManager(obj,iC,lG,bS,mR,iCR)
            s.nSubdomains     = obj.nSubdomains;
            s.interfaceConnec = iC;
            s.interfaceConnecReshaped = iCR;
            s.locGlobConnec   = lG;
            s.nBoundaryNodes  = bS{1}.mesh.nnodes;
            s.nReferenceNodes = mR.nnodes;
            s.nNodes          = obj.mesh.nnodes;
            s.nDimf           = obj.mesh.ndim;
            d = DomainDecompositionDofManager(s);
        end

        function c = createComplianceFromConstiutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstitutiveTensor(s);
        end

        function createComplianceRadius(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.EIFEMprecontitioner;
            obj.compliance = ComplianceFunctionalRadius(s);
        end

        function createCompliance(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filter;
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function uMesh = createBaseDomain(obj)
%             levelSet         = -ones(obj.mesh.nnodes,1);
%             s.backgroundMesh = obj.mesh;
%             s.boundaryMesh   = obj.mesh.createBoundaryMesh();
%             uMesh = UnfittedMesh(s);
%             uMesh.compute(levelSet);
              uMesh = obj.designVariable.fun;
        end
        
        function createVolumeConstraint(obj)
            s.mesh   = obj.coarseMesh;
            s.filter = obj.filter;
            s.test = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = obj.volumeTarget;
            s.uMesh = obj.createBaseDomain();
            s.geomType = 'Circle';
            v = VolumeConstraintRadius(s);
            obj.volume = v;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            test   = LagrangianFunction.create(obj.coarseMesh, 1, 'P0');
            trial  = LagrangianFunction.create(obj.coarseMesh, 1, 'P0');
            f = @(u,v) DP(v,u);
            M = IntegrateLHS(f,test,trial,obj.coarseMesh,'Domain',2);
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createDualVariable(obj)
            s.nConstraints   = 1;
            l                = DualVariable(s);
            obj.dualVariable = l;
        end

         function createPrimalUpdater(obj)
            s.ub     = 0.96;
            s.lb     = 1e-6; % fro lattice 0.01;
            s.tauMax = 1000;
            s.tau    = [];
            obj.primalUpdater = ProjectedGradient(s);
        end

        function createOptimizer(obj)
%             s.monitoring     = true;
%             s.cost           = obj.cost;
%             s.constraint     = obj.constraint;
%             s.designVariable = obj.designVariable;
%             s.dualVariable   = obj.dualVariable;
%             s.maxIter        = 3000;
%             s.tolerance      = 1e-8;
%             s.constraintCase = {'EQUALITY'};
%             s.ub             = 0.8;
%             s.lb             = 1e-6;
%             s.volumeTarget   = obj.volumeTarget;
%             s.primal         = 'PROJECTED GRADIENT';
%             opt              = OptimizerMMA(s);
%             opt.solveProblem();
%             obj.optimizer = opt;
           
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 2000; %150
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.primal         = 'PROJECTED GRADIENT';
            s.etaNorm        = 10;%0.05
            s.gJFlowRatio    = 10; %3
            s.primalUpdater  = obj.primalUpdater;
%             s.etaMaxMin      = 0.05;
%             s.etaMax      = 200;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj)
            xMin    = min(obj.mesh.coord(:,1));
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            yMin    = min(obj.mesh.coord(:,2)); 
            Ly      = yMax-yMin;
            isDir   = @(coor)  abs(coor(:,1)-xMin) < 1e-12;
            isForce = @(x)  abs(x(1,:,:)-xMax) < 1e-12 & x(2,:,:)>=yMin+0.3*Ly & x(2,:,:)<=yMin+0.7*Ly;

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            [bMesh, ~]  = obj.mesh.createSingleBoundaryMesh();
            sPL{1}.domain = isForce;
            sPL{1}.fun    = ConstantFunction.create([0,-1],bMesh);

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = TractionLoad(obj.mesh,sPL{i},'FUNCTION');
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
            obj.sDir = sDir;
        end
    end
end
