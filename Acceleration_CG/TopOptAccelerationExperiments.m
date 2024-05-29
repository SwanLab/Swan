classdef TopOptAccelerationExperiments < handle

    properties (Access = public)
        mesh
        filter
        designVariable
        materialInterpolator
        physicalProblem
        solverTol
        compliance
        volume
        cost
        constraint
        dualVariable
        optimizer
        momentum
    end

    methods (Access = public)

        function obj = TopOptAccelerationExperiments()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createVolumeConstraint();
            obj.createSolverTolerance();
            obj.createElasticProblem();
            obj.createComplianceFromConstitutive();
            obj.createCompliance();
            obj.createCost();
            obj.createConstraint();
            obj.createDualVariable();
            obj.createMomentum();
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            %UnitMesh better
            % 2D cantilever
            x1      = linspace(0,2,100);
            x2      = linspace(0,1,50);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh.create(s);

            % 2D arch

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
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = f.project('P1');            
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            m = Material.create(s);
        end

        function createSolverTolerance(obj)
            s.solver      = 'CONJUGATE GRADIENT';
            s.tolMax      = 2e-1;
            s.tolMin      = 1e-5;
            obj.solverTol = ConjugateGradientToleranceCalculator(s);           
        end

        function createElasticProblem(obj)
            s.mesh     = obj.mesh;
            s.scale    = 'MACRO';
            s.material = obj.createMaterial();
            s.dim      = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType  = 'LINEAR';
            s.solverType         = 'REDUCED';
            s.solverMode         = 'DISP';
            s.solverCase         = 'CONJUGATE GRADIENT';%'DIRECT';%
            s.solverTol          = obj.solverTol;
            s.matrixFree         = false;
            p.maxIters           = 5e3;
            p.displayInfo        = true;
            s.solverParams       = p;
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function c = createComplianceFromConstitutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstiutiveTensor(s);
        end

        function createCompliance(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filter;
            s.complainceFromConstitutive  = obj.createComplianceFromConstitutive();
            s.material                    = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.2;
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
            s.test  = LagrangianFunction.create(obj.mesh,1,'P1');
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh  = obj.mesh;
            s.type  = 'MassMatrix';
            LHS = LHSintegrator.create(s);
            M = LHS.compute;     
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

        function createMomentum(obj)
            s.momentumCase = 'Nesterov';
            s.betaStrategy = 'Constant';
            s.momentumVal  = 0;
            s.x0           = obj.designVariable.fun.fValues;
            obj.momentum   = Momentum(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 1000;
            s.tolerance      = 5e-2;
            s.constraintCase = {'EQUALITY'};
            s.ub             = 1;
            s.lb             = 0;
            s.volumeTarget   = 0.1;
            s.primal         = 'PROJECTED GRADIENT';
            s.solverTol      = obj.solverTol;
            s.constantTau    = true;
            s.tauValue       = 5e-3;%5e-2;%1e-2;%
            s.momentum       = obj.momentum;
            s.etaNorm        = 100;
            s.gJFlowRatio    = 100;
            s.etaMax         = 1.5e3;
            % opt = OptimizerAugmentedLagrangian(s);
            opt = OptimizerNullSpace(s);
            % opt = OptimizerMMA(s);
            tStart = tic;
            opt.solveProblem();
            toc(tStart)
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));

            % cantilever
            % isDir   = @(coor)  abs(coor(:,1))==0;
            % isForce = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.3*yMax & abs(coor(:,2))<=0.7*yMax);

            % arch
            isDir   = @(coor) abs(coor(:,2)) == 0 & (abs(coor(:,1))<=0.2*xMax | abs(coor(:,1))>=0.8*xMax);
            isForce = @(coor) abs(coor(:,2)) == 0 & abs(coor(:,1))>=0.4*xMax & abs(coor(:,1))<=0.6*xMax; 

            % bridge
            % isDir   = @(coor) abs(coor(:,2)) == 0;
            % isForce = @(coor) abs(coor(:,2)) == yMax;



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
                pl = PointLoad(obj.mesh, sPL{i});
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
        end
    end
end
