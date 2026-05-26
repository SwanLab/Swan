classdef TutorialFirst < handle

    properties (Access = public)
        designVariable  % <-- público para exportar b* para Fase 2
    end

    properties (Access = private)
        mesh
        physicalProblem
        materialMicro
        cost
        optimizer
        filterRegularization
    end

    methods (Access = public)

        function obj = TutorialFirst()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilterRegularization();
            obj.createMaterial();           % <-- antes do elastic problem
            obj.createElasticProblem();
            obj.createBaseDomain();
            obj.createCost();
            obj.createOptimizer();          % <-- sem constraint
        end
        function hist = getCostHistory(obj)
            hist = obj.optimizer.costHistory;
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            obj.mesh = TriangleMesh(2, 1, 100, 60);
        end

        function createDesignVariable(obj)
            % Só b — rho=1 fixo (material cheio)
            s.fHandle = @(x) zeros(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);

            sD.fun      = aFun.project('P1');
            sD.mesh     = obj.mesh;
            sD.type     = 'Density';
            sD.plotting = true;
            obj.designVariable = DesignVariable.create(sD);
        end

        function createFilterRegularization(obj)
            eOverhmin = 6;
            s.filterType   = 'PDE';
            s.mesh         = obj.mesh;
            s.boundaryType = 'Neumann';
            s.metric       = 'Isotropy';
            s.trial        = LagrangianFunction.create(obj.mesh, 1, 'P1');
            f = Filter.create(s);
            epsilon = eOverhmin * obj.mesh.computeMeanCellSize();
            f.updateEpsilon(epsilon);
            obj.filterRegularization = f;
        end

        function createMaterial(obj)
            % Só C^h(b) — sem SIMP, rho=1
            s.type     = 'HomogenizedMicrostructure';
            s.mesh     = obj.mesh;
            s.young    = 1.0;
            s.fileName = 'HomogenizationLattice4';
            s.density  = obj.designVariable;
            obj.materialMicro = MaterialFactory.create(s);
        end

        function createElasticProblem(obj)
            s.mesh               = obj.mesh;
            s.scale              = 'MACRO';
            s.material           = obj.materialMicro;
            s.dim                = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType  = 'LINEAR';
            s.solverType         = 'REDUCED';
            s.solverMode         = 'DISP';
            s.solverCase         = DirectSolver();
            obj.physicalProblem  = ElasticProblem(s);
        end

        function uMesh = createBaseDomain(obj)
            levelSet         = -ones(obj.mesh.nnodes, 1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh            = UnfittedMesh(s);
            uMesh.compute(levelSet);
        end

        function c = createComplianceFromConstitutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstitutiveTensor(s);
        end

        function c = createCompliance(obj)
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filterRegularization;
            s.complainceFromConstitutive = obj.createComplianceFromConstitutive();
            s.material                   = obj.materialMicro;
            c = ComplianceFunctional(s);  % <-- funcional simples, sem densidade
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.createCompliance();
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            test  = LagrangianFunction.create(obj.mesh, 1, 'P1');
            trial = LagrangianFunction.create(obj.mesh, 1, 'P1');
            M = IntegrateLHS(@(u,v) DP(v,u), test, trial, obj.mesh, 'Domain');
        end

        function createOptimizer(obj)
            s.cost           = obj.cost;
            s.designVariable = obj.designVariable;
            s.monitoring     = true;
            s.lb             = -0.6;
            s.ub             =  0.6;
            s.maxIter        = 1000;
            opt              = OptimizerProjectedGradient(s);
            opt.solveProblem();
            obj.optimizer    = opt;
        end

        function bc = createBoundaryConditions(obj)
            xMax = max(obj.mesh.coord(:,1));
            yMax = max(obj.mesh.coord(:,2));
            isDir   = @(coor) abs(coor(:,1)) == 0;
            isForce = @(coor) abs(coor(:,1)) == xMax & ...
                              abs(coor(:,2)) >= 0.4*yMax & ...
                              abs(coor(:,2)) <= 0.6*yMax;

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1, 2];
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