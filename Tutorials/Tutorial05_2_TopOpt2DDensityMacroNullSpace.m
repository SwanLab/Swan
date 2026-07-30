classdef Tutorial05_2_TopOpt2DDensityMacroNullSpace < handle

    properties (Access = private)
        mesh
        filter
        designVariable
        C
        dC
        physicalProblem
        compliance
        volume
        cost
        constraint
        primalUpdater
        optimizer
        E0; nu0; E1; nu1
    end

    methods (Access = public)

        function obj = Tutorial05_2_TopOpt2DDensityMacroNullSpace()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterial();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
            obj.E0  = 1e-3;
            obj.nu0 = 1/3;
            obj.E1  = 1;
            obj.nu1 = 1/3;
        end

        function createMesh(obj)
            obj.mesh = TriangleMesh(2,1,100,50);
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);

            sD.fun      = aFun.project('P1');
            sD.mesh     = obj.mesh;
            sD.type     = 'Density';
            sD.plotting = true;
            dens        = DesignVariable.create(sD);
            obj.designVariable = dens;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

        function createMaterial(obj)
            N = obj.mesh.ndim;
            muA    = obj.computeMu(obj.E0,obj.nu0);
            kappaA = obj.computeKappa(obj.E0,obj.nu0,N);
            muB    = obj.computeMu(obj.E1,obj.nu1);
            kappaB = obj.computeKappa(obj.E1,obj.nu1,N);

            mu    = @(rho) SimpAllInterpolator.computeMu(muA,muB,kappaA,kappaB,rho,N); mu = @(rho) Expand(mu(rho),4);
            kappa  = @(rho) SimpAllInterpolator.computeKappa(muA,muB,kappaA,kappaB,rho,N); kappa = @(rho) Expand(kappa(rho),4);
            lambda = @(rho) kappa(rho) - (2/N)*mu(rho);
            I      = ConstantFunction.create(eye4D(N),obj.mesh);
            IxI    = ConstantFunction.create(kronEye(N),obj.mesh);
            obj.C  = @(rho) 2*mu(rho).*I + lambda(rho).*IxI;

            dmu     = @(rho) SimpAllInterpolator.computeMuDerivative(muA,muB,kappaA,kappaB,rho,N); dmu = @(rho) Expand(dmu(rho),4);
            dkappa  = @(rho) SimpAllInterpolator.computeKappaDerivative(muA,muB,kappaA,kappaB,rho,N); dkappa = @(rho) Expand(dkappa(rho),4);
            dlambda = @(rho) dkappa(rho) - (2/N)*dmu(rho);
            obj.dC  = @(rho) 2*dmu(rho).*I + dlambda(rho).*IxI;
        end

        function mu = computeMu(obj,E,nu)
            mu = E./(2*(1+nu));
        end

        function kappa = computeKappa(obj,E,nu,N)
            kappa = E./(N*(1-(N-1)*nu));
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = [];
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = DirectSolver();
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
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
            s.C                           = obj.C;
            s.dC                          = obj.dC;
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
            test  = LagrangianFunction.create(obj.mesh,1,'P1');
            trial = LagrangianFunction.create(obj.mesh,1,'P1');
            M = IntegrateLHS(@(u,v) DP(v,u),test,trial,obj.mesh,'Domain');
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createPrimalUpdater(obj)
            s.ub     = 1;
            s.lb     = 0;
            s.tauMax = 1000;
            s.tau    = [];
            obj.primalUpdater = ProjectedGradient(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 3;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.primal         = 'PROJECTED GRADIENT';
            s.delta          = 0.01;
            s.etaStar        = 2;
            s.primalUpdater  = obj.primalUpdater;
            s.gif            = false;
            s.gifName        = [];
            s.printing       = false;
            s.printName      = [];
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.4*yMax & abs(coor(:,2))<=0.6*yMax);

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = -1;

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = TractionLoad(obj.mesh, sPL{i}, 'DIRAC');
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
        end
    end
end