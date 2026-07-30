classdef Tutorial05_1_TopOpt2DDensityMacroMMA < handle

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
        E0; nu0; E1; nu1
    end

    methods (Access = public)

        function obj = Tutorial05_1_TopOpt2DDensityMacroMMA()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialProperties();
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
        end

        function createMesh(obj)
            obj.mesh = TriangleMesh(2,1,100,50);
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);
            s.fun     = aFun.project('P1');
            s.mesh    = obj.mesh;
            s.type = 'Density';
            s.plotting = true;
            dens    = DesignVariable.create(s);
            obj.designVariable = dens;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

        function createMaterialProperties(obj)
            obj.E0 = 1e-3;
            obj.nu0 = 1/3;
            obj.E1 = 1;
            obj.nu1 = 1/3;
        end

        function [C,dC] = createMaterial(obj)
            N = obj.mesh.ndim;
            muA    = obj.computeMu(obj.E0,obj.nu0);
            kappaA = obj.computeKappa(obj.E0,obj.nu0,N);
            muB    = obj.computeMu(obj.E1,obj.nu1);
            kappaB = obj.computeKappa(obj.E1,obj.nu1,N);

            mu    = @(rho) SimpP3Interpolator.computeMu(muA,muB,rho); mu = @(rho) Expand(mu(rho),4);
            kappa  = @(rho) SimpP3Interpolator.computeKappa(kappaA,kappaB,rho); kappa = @(rho) Expand(kappa(rho),4);
            lambda = @(rho) kappa(rho) - (2/N)*mu(rho);
            I      = ConstantFunction.create(eye4D(N),obj.mesh);
            IxI    = ConstantFunction.create(kronEye(N),obj.mesh);
            C  = @(rho) 2*mu(rho).*I + lambda(rho).*IxI;

            dmu     = @(rho) SimpP3Interpolator.computeMuDerivative(muA,muB,rho); dmu = @(rho) Expand(dmu(rho),4);
            dkappa  = @(rho) SimpP3Interpolator.computeKappaDerivative(kappaA,kappaB,rho); dkappa = @(rho) Expand(dkappa(rho),4);
            dlambda = @(rho) dkappa(rho) - (2/N)*dmu(rho);
            dC  = @(rho) 2*dmu(rho).*I + dlambda(rho).*IxI;
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
            [C,dC] = obj.createMaterial();
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filter;
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.C                           = C;
            s.dC                          = dC;
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
            M = diag(sum(M,1));
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
            s.gif            = true;
            s.gifName        = 'Tutorial05_1';
            s.printing       = true;
            s.printName      = 'Tutorial05_1';
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
