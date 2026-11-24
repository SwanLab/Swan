classdef Tutorial_TopOpt2DLatticeTurnTMacroNNNullSpace < handle

    properties (Access = private)
        mesh
        filters
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
        compliance0
    end

    methods (Access = public)

        function obj = Tutorial_TopOpt2DLatticeTurnTMacroNNNullSpace()

            obj.init();
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilters();
            obj.createMaterial();
            obj.createElasticProblem();
            obj.createCompliance();
            obj.computeInitialCompliance();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();

        end

        function J = getComplianceHistory(obj)

            J = obj.compliance0*obj.cost.getIterationHistory();

        end

        function saveDesign(obj,fileName)

            fields = obj.designVariable.fields;
            mU     = latticeXiSettings().hole.mU;

            rho = obj.filters{1}.compute(fields{1},2).fValues;
            a   = obj.filters{2}.compute(fields{2},2).fValues;
            t   = obj.filters{3}.compute(fields{3},2).fValues;

            [xiMin,xiMax] = XiLatticeDatasetGenerator.computeXiBounds( ...
                sqrt(1-rho),mU);

            design.coord  = obj.mesh.coord;
            design.connec = obj.mesh.connec;
            design.rho    = rho;
            design.a      = a;
            design.t      = t;
            % design.xi     = obj.filters{3}.compute(fields{3},2).fValues;
            design.xi     = xiMin + t.*(xiMax-xiMin);
            design.mU     = mU;

            save(fileName,'-struct','design');

        end

        function plotFields(obj)

            fields = obj.designVariable.fields;
            limits = obj.createFieldLimits();

            plotLatticeFields(fields,{'\rho','a','t'},limits);

            rho = squeeze(fields{1}.project('P0').fValues);
            t   = squeeze(fields{3}.project('P0').fValues);
            mU  = latticeXiSettings().hole.mU;

            [xiMin,xiMax] = XiLatticeDatasetGenerator.computeXiBounds( ...
                sqrt(1-rho),mU);

            xi = xiMin + t.*(xiMax-xiMin);

            plotFieldPatch(obj.mesh.connec,obj.mesh.coord,xi,'\xi',[0,pi/2]);

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

            sRho.fHandle = @(x) 0.99*ones(size(x(1,:,:)));
            sRho.ndimf   = 1;
            sRho.mesh    = obj.mesh;

            rhoAnalytical = AnalyticalFunction(sRho);
            rho = rhoAnalytical.project('P1');

            sA.fHandle = @(x) zeros(size(x(1,:,:)));
            sA.ndimf   = 1;
            sA.mesh    = obj.mesh;

            aAnalytical = AnalyticalFunction(sA);
            a = aAnalytical.project('P1');

            sT.fHandle = @(x) 0.5*ones(size(x(1,:,:)));
            sT.ndimf   = 1;
            sT.mesh    = obj.mesh;

            tAnalytical = AnalyticalFunction(sT);
            t = tAnalytical.project('P1');

            sD.fields   = {rho,a,t};
            sD.type     = 'MultiField';
            sD.plotting = true;

            obj.designVariable = DesignVariable.create(sD);

        end

        function createFilters(obj)

            h = obj.mesh.computeMeanCellSize();

            sRho.filterType = 'LUMP';
            sRho.mesh       = obj.mesh;
            sRho.trial      = LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            filterRho = Filter.create(sRho);

            sA.filterType   = 'PDE';
            sA.mesh         = obj.mesh;
            sA.boundaryType = 'Neumann';
            sA.metric       = 'Isotropy';
            sA.trial        = LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            filterA = Filter.create(sA);
            filterA.updateEpsilon(3*h);

            % sT.filterType   = 'PDE';
            % sT.mesh         = obj.mesh;
            % sT.boundaryType = 'Neumann';
            % sT.metric       = 'Isotropy';
            % sT.trial        = LagrangianFunction.create( ...
            %     obj.mesh,1,'P1');
            %
            % filterT = Filter.create(sT);
            % filterT.updateEpsilon(3*h);

            sT.filterType = 'LUMP';
            sT.mesh       = obj.mesh;
            sT.trial      = LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            filterT = Filter.create(sT);

            obj.filters = {filterRho,filterA,filterT};

        end

        function createMaterial(obj)

            % s.fileName      = 'latticeXiTNN_Adam_64x3_cholesky';
            s.fileName      = 'latticeSquareTurnTNN_Adam_64x3_cholesky';
            s.mesh          = obj.mesh;
            s.young         = 1;
            s.parameterNames = {'rho','a','t'};

            hm = LatticeHomogenizedMaterialsReader(s);

            obj.C = @(designVariables) ...
                hm.obtainTensor(designVariables);

            obj.dC = @(designVariables) ...
                hm.obtainTensorDerivatives(designVariables);

        end

        function createElasticProblem(obj)

            s.mesh               = obj.mesh;
            s.scale              = 'MACRO';
            s.material           = [];
            s.dim                = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType  = 'LINEAR';
            s.solverType         = 'REDUCED';
            s.solverMode         = 'DISP';
            s.solverCase         = DirectSolver();

            obj.physicalProblem = ElasticProblem(s);

        end

        function c = createComplianceFromConstitutive(obj)

            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;

            c = ComplianceFromConstitutiveTensor(s);

        end

        function createCompliance(obj)

            s.mesh                       = obj.mesh;
            s.filters                    = obj.filters;
            s.complainceFromConstitutive = ...
                obj.createComplianceFromConstitutive();
            s.C                          = obj.C;
            s.dC                         = obj.dC;

            obj.compliance = ComplianceFunctional(s);

        end

        function computeInitialCompliance(obj)

            s.mesh                       = obj.mesh;
            s.filters                    = obj.filters;
            s.complainceFromConstitutive = ...
                obj.createComplianceFromConstitutive();
            s.C                          = obj.C;
            s.dC                         = obj.dC;
            s.value0                     = 1;

            physicalCompliance = ComplianceFunctional(s);

            obj.compliance0 = physicalCompliance.computeFunctionAndGradient( ...
                obj.designVariable);

        end

        function uMesh = createBaseDomain(obj)

            levelSet = -ones(obj.mesh.nnodes,1);

            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();

            uMesh = UnfittedMesh(s);

            uMesh.compute(levelSet);

        end

        function createVolumeConstraint(obj)

            s.mesh             = obj.mesh;
            s.test             = LagrangianFunction.create( ...
                obj.mesh,1,'P1');
            s.volumeTarget     = 0.4;
            s.volumeFieldIndex = 1;
            s.uMesh            = obj.createBaseDomain();

            obj.volume = VolumeConstraint(s);

        end

        function createCost(obj)

            s.shapeFunctions{1} = obj.compliance;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();

            obj.cost = RecordingCost(Cost(s));

        end

        function M = createMassMatrix(obj)

            test = LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            trial = LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            MSingle = IntegrateLHS( ...
                @(u,v) DP(v,u), ...
                test, ...
                trial, ...
                obj.mesh, ...
                'Domain');

            MSingle = diag(sum(MSingle,1));

            nFields = obj.designVariable.obtainNumberVariables();

            M = kron(speye(nFields),MSingle);

        end

        function createConstraint(obj)

            s.shapeFunctions{1} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();

            obj.constraint = Constraint(s);

        end

        function limits = createFieldLimits(obj)

            rhoMin = 1e-8;
            rhoMax = 0.99;
            % aMax   = 0.8;
            aMax   = 1;

            limits = {[rhoMin,rhoMax],[-aMax,aMax],[0,1]};

        end

        function createPrimalUpdater(obj)

            fields = obj.designVariable.obtainDomainFunction();
            limits = obj.createFieldLimits();

            s.lb = [];
            s.ub = [];

            for iField = 1:numel(fields)

                n    = numel(fields{iField}.fValues);
                s.lb = [s.lb; limits{iField}(1)*ones(n,1)];
                s.ub = [s.ub; limits{iField}(2)*ones(n,1)];

            end

            s.tauMax = 1000;
            s.tau    = [];

            obj.primalUpdater = ProjectedGradient(s);

        end

        function createOptimizer(obj)

            s.monitoring      = true;
            s.cost            = obj.cost;
            s.constraint      = obj.constraint;
            s.designVariable  = obj.designVariable;
            s.maxIter         = 2000;
            s.tolerance       = 1e-8;
            s.constraintCase  = {'EQUALITY'};
            s.primal          = 'PROJECTED GRADIENT';
            s.delta           = 0.01;
            s.etaStar         = 2;
            s.primalUpdater   = obj.primalUpdater;
            s.gif             = false;
            s.gifName         = [];
            s.printing        = false;
            s.printName       = [];

            opt = OptimizerNullSpace(s);

            opt.solveProblem();

            obj.optimizer = opt;

        end

        function bc = createBoundaryConditions(obj)

            xMax = max(obj.mesh.coord(:,1));
            yMax = max(obj.mesh.coord(:,2));

            isDir = @(coor) ...
                abs(coor(:,1)) == 0;

            isForce = @(coor) ...
                (abs(coor(:,1)) == xMax & ...
                abs(coor(:,2)) >= 0.4*yMax & ...
                abs(coor(:,2)) <= 0.6*yMax);

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = -1;

            dirichletFun = [];

            for i = 1:numel(sDir)

                dir = DirichletCondition( ...
                    obj.mesh, ...
                    sDir{i});

                dirichletFun = [dirichletFun,dir];

            end

            s.dirichletFun = dirichletFun;

            pointloadFun = [];

            for i = 1:numel(sPL)

                pl = TractionLoad( ...
                    obj.mesh, ...
                    sPL{i}, ...
                    'DIRAC');

                pointloadFun = [pointloadFun,pl];

            end

            s.pointloadFun = pointloadFun;
            s.periodicFun  = [];
            s.mesh         = obj.mesh;

            bc = BoundaryConditions(s);

        end

    end

end