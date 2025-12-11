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
        fileNameEIFEM
    end

    methods (Access = public)

        function obj = OptimizationEIFEMTutorial()
            obj.init()
            obj.createMesh();
            obj.createCoarseMesh();
            obj.createDesignVariable();
%             obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createDualVariable();
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
            obj.nSubdomains = [5,3];
            obj.r = 1e-6*ones(obj.nSubdomains)'; 
            obj.r= (0.6 - 0.2) * rand(obj.nSubdomains(2),obj.nSubdomains(1)) + 0.1;
            obj.xmax=1; obj.xmin=-1; obj.ymax = 1; obj.ymin=-1; 
            obj.Nr = 7; obj.Ntheta = 14; 
            obj.x0 = 0; obj.y0=0;
            obj.tolSameNode = 1e-10;
            obj.fileNameEIFEM = './EPFL/parametrizedEIFEMLagrange40.mat';
        end

        function createMesh(obj)
            mSbd = obj.createSubDomainMeshes();
            bS = mSbd{1,1}.createBoundaryMesh();
            [mD,mSb,iC,lG,iCR,discMesh] = obj.createMeshDomainJoiner(mSbd);
            obj.mesh = mD;
            obj.subdomainMeshes = mSb;
        end

         function  mSbd = createSubDomainMeshes(obj)
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            Lx = obj.xmax-obj.xmin;
            Ly = obj.ymax-obj.ymin;
            for jDom = 1:nY
                for iDom = 1:nX
                    refMesh = mesh_rectangle_via_triangles(obj.r(jDom,iDom),obj.xmax,obj.xmin,obj.ymax,obj.ymin,obj.Nr,obj.Ntheta,obj.x0,obj.y0);
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
            %                     connec = [1 2 3 4];
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
            s.plotting = false;
            radius    = DesignVariable.create(s);
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

        function m = createMaterial(obj)
            x = obj.designVariable.fun;           
            s.type                 = 'DensityBased';
            s.density              = x;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            s.mesh                 = obj.mesh;
            m = Material.create(s);
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function Meifem = createEIFEMPreconditioner(obj,dir,iC,lG,bS,iCR,dMesh,mSbd)
            mR            = obj.referenceMesh;
            s.RVE         = TrainedRVE(obj.fileNameEIFEM);
            s.mesh        = obj.createCoarseMesh(mR);
            s.DirCond     = dir;
            s.nSubdomains = obj.nSubdomains;
            s.mu          = obj.r;
            s.meshRef     = dMesh;
            eifem         = EIFEMnonPeriodic(s);

            ss.ddDofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,mR,iCR);
            ss.EIFEMsolver = eifem;
            ss.bcApplier = obj.bcApplier;
            ss.dMesh     = dMesh;
            ss.type = 'EIFEM';
            eP = Preconditioner.create(ss);
            Meifem = @(r) eP.apply(r);
        end

        function c = createComplianceFromConstiutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstitutiveTensor(s);
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
            levelSet         = -ones(obj.mesh.nnodes,1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
            uMesh.compute(levelSet);
        end
        
        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.test = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;
            s.uMesh = obj.createBaseDomain();
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            test   = LagrangianFunction.create(obj.mesh, 1, 'P1');
            trial  = LagrangianFunction.create(obj.mesh, 1, 'P1');
            f = @(u,v) DP(v,u);
            M = IntegrateLHS(f,test,trial,obj.mesh,'Domain',2);
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

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 3;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.ub             = 1;
            s.lb             = 0;
            s.volumeTarget   = 0.4;
            s.primal         = 'PROJECTED GRADIENT';
            opt              = OptimizerMMA(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(x)  x(1,:,:)==xMax & x(2,:,:)>=0.3*yMax & x(2,:,:)<=0.7*yMax;

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
        end
    end
end
