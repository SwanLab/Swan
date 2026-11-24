classdef Tutorial_TopOpt2DLatticeSidesMacroNNNullSpace < handle

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
        initialDesignFile
        initialVolumeFraction
    end

    methods (Access = public)

        function obj = Tutorial_TopOpt2DLatticeSidesMacroNNNullSpace(runOptimization,initialDesignFile,initialVolumeFraction)

            if nargin < 1 || isempty(runOptimization)
                runOptimization = true;
            end
            if nargin < 2
                initialDesignFile = '';
            end
            if nargin < 3
                initialVolumeFraction = [];
            end

            obj.initialDesignFile     = initialDesignFile;
            obj.initialVolumeFraction = initialVolumeFraction;

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

            if runOptimization
                obj.createOptimizer();
            end

        end

        function J = getComplianceHistory(obj)

            J = obj.compliance0*obj.cost.getIterationHistory();

        end

        function saveDesign(obj,fileName)

            [m1,m2,a] = obj.computeFilteredFields();

            design.coord  = obj.mesh.coord;
            design.connec = obj.mesh.connec;
            design.m1     = m1;
            design.m2     = m2;
            design.a      = a;
            design.rho    = 1 - m1.*m2;
            design.xi     = atan2(m1,m2);
            design.mU     = latticeXiSettings().hole.mU;

            save(fileName,'-struct','design');

        end

        function plotFields(obj)

            fields = obj.designVariable.fields;
            limits = obj.createFieldLimits();

            plotLatticeFields(fields,{'m_1','m_2','a'},limits);

            m1  = squeeze(fields{1}.project('P0').fValues);
            m2  = squeeze(fields{2}.project('P0').fValues);
            rho = 1 - m1.*m2;
            xi  = atan2(m1,m2);

            plotFieldPatch(obj.mesh.connec,obj.mesh.coord,rho,'\rho',[0,1]);
            plotFieldPatch(obj.mesh.connec,obj.mesh.coord,xi,'\xi',[0,pi/2]);

        end

        function err = checkGradients(obj,h,seed,direction)

            if nargin < 2 || isempty(h)
                h = 1e-5;
            end
            if nargin < 3 || isempty(seed)
                seed = 1;
            end
            if nargin < 4 || isempty(direction)
                direction = 'smooth';
            end

            rng(seed);

            x0     = obj.designVariable.fun.fValues;
            fields = obj.designVariable.fields;
            sizes  = cellfun(@(f) numel(f.fValues),fields);
            offset = [0,cumsum(sizes(:)')];

            obj.cost.computeFunctionAndGradient(obj.designVariable);
            obj.constraint.computeFunctionAndGradient(obj.designVariable);
            gJ = obj.cost.gradient;
            gG = obj.constraint.gradient;

            err = zeros(numel(fields),2);

            for i = 1:numel(fields)

                d = zeros(size(x0));
                d(offset(i)+1:offset(i+1)) = obj.createTestDirection(direction);

                obj.designVariable.update(x0 + h*d);
                obj.cost.computeFunctionAndGradient(obj.designVariable);
                obj.constraint.computeFunctionAndGradient(obj.designVariable);
                Jp = obj.cost.value;
                Gp = obj.constraint.value;

                obj.designVariable.update(x0 - h*d);
                obj.cost.computeFunctionAndGradient(obj.designVariable);
                obj.constraint.computeFunctionAndGradient(obj.designVariable);
                Jm = obj.cost.value;
                Gm = obj.constraint.value;

                fdJ = (Jp - Jm)/(2*h);
                fdG = (Gp - Gm)/(2*h);
                anJ = gJ'*d;
                anG = gG'*d;

                err(i,:) = [abs(fdJ-anJ)/max(abs(anJ),eps), ...
                    abs(fdG-anG)/max(abs(anG),eps)];

            end

            obj.designVariable.update(x0);
            obj.cost.computeFunctionAndGradient(obj.designVariable);
            obj.constraint.computeFunctionAndGradient(obj.designVariable);

            disp(table(err(:,1),err(:,2),'RowNames',{'m1','m2','a'}, ...
                'VariableNames',{'cost','volume'}))

        end

        function d = createTestDirection(obj,direction)

            x  = obj.mesh.coord(:,1);
            y  = obj.mesh.coord(:,2);
            Lx = max(x) - min(x);
            Ly = max(y) - min(y);

            switch direction
                case 'random'
                    d = rand(obj.mesh.nnodes,1) - 0.5;
                case 'smooth'
                    d = zeros(obj.mesh.nnodes,1);
                    for k = 1:3
                        c = rand(1,4) - 0.5;
                        d = d + c(1)*sin(k*pi*x/Lx).*sin(k*pi*y/Ly) ...
                            + c(2)*cos(k*pi*x/Lx).*cos(k*pi*y/Ly) ...
                            + c(3)*sin(k*pi*x/Lx).*cos(k*pi*y/Ly) ...
                            + c(4)*cos(k*pi*x/Lx).*sin(k*pi*y/Ly);
                    end
                otherwise
                    error('checkGradients:UnknownDirection','Unknown direction %s.',direction);
            end

        end
    end

    methods (Access = private)

        function init(obj)

            close all;

        end

        function createMesh(obj)

            obj.mesh = TriangleMesh(2,1,100,50);

        end

        function [mL,mU,aMax] = createDesignParameters(obj)

            mL   = 0.1;
            mU   = latticeXiSettings().hole.mU;
            aMax = 0.8;

        end

        function createDesignVariable(obj)

            if isempty(obj.initialDesignFile)
                [m1,m2,a] = obj.createUniformInitialFields();
            else
                [m1,m2,a] = obj.loadInitialFields();
            end

            sD.fields   = {m1,m2,a};
            sD.type     = 'MultiField';
            sD.plotting = true;

            obj.designVariable = DesignVariable.create(sD);

        end

        function [m1,m2,a] = createUniformInitialFields(obj)

            [mL,mU,~] = obj.createDesignParameters();

            % m0 = mL;
            if isempty(obj.initialVolumeFraction)
                m0 = mL;
            else
                m0 = min(mU,max(mL,sqrt(1 - obj.initialVolumeFraction)));
            end

            sM.fHandle = @(x) m0*ones(size(x(1,:,:)));
            sM.ndimf   = 1;
            sM.mesh    = obj.mesh;

            m1 = AnalyticalFunction(sM).project('P1');
            m2 = AnalyticalFunction(sM).project('P1');

            sA.fHandle = @(x) zeros(size(x(1,:,:)));
            sA.ndimf   = 1;
            sA.mesh    = obj.mesh;

            a = AnalyticalFunction(sA).project('P1');

        end

        function [m1,m2,a] = loadInitialFields(obj)

            [mL,mU,aMax] = obj.createDesignParameters();

            design = load(obj.initialDesignFile);

            if ~isequal(size(design.coord),size(obj.mesh.coord)) || ...
                    max(abs(design.coord(:) - obj.mesh.coord(:))) > 1e-12
                error('Tutorial_TopOpt2DLatticeSidesMacroNNNullSpace:MeshMismatch', ...
                    'The mesh of %s does not match the tutorial mesh.',obj.initialDesignFile);
            end

            holeArea = 1 - design.rho(:);
            tanXi    = tan(design.xi(:));

            m1Values = min(mU,max(mL,sqrt(holeArea.*tanXi)));
            m2Values = min(mU,max(mL,sqrt(holeArea./tanXi)));
            aValues  = min(aMax,max(-aMax,design.a(:)));

            m1 = LagrangianFunction.create(obj.mesh,1,'P1');
            m2 = LagrangianFunction.create(obj.mesh,1,'P1');
            a  = LagrangianFunction.create(obj.mesh,1,'P1');

            m1.setFValues(m1Values);
            m2.setFValues(m2Values);
            a.setFValues(aValues);

        end

        function createFilters(obj)

            h     = obj.mesh.computeMeanCellSize();
            types = {'LUMP','LUMP','PDE'};

            obj.filters = cell(1,numel(types));

            for i = 1:numel(types)

                s       = struct();
                s.mesh  = obj.mesh;
                s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
                s.filterType = types{i};

                if strcmp(types{i},'PDE')
                    s.boundaryType = 'Neumann';
                    s.metric       = 'Isotropy';
                end

                f = Filter.create(s);

                if strcmp(types{i},'PDE')
                    f.updateEpsilon(3*h);
                end

                obj.filters{i} = f;

            end

        end

        function [m1,m2,a] = computeFilteredFields(obj)

            fields = obj.designVariable.fields;

            m1 = obj.filters{1}.compute(fields{1},2).fValues;
            m2 = obj.filters{2}.compute(fields{2},2).fValues;
            a  = obj.filters{3}.compute(fields{3},2).fValues;

        end

        function createMaterial(obj)

            % s.fileName       = 'latticeSidesNN_Adam_64x3_cholesky';
            s.fileName       = 'latticeSidesNN_native_Adam_64x3_cholesky';
            s.mesh           = obj.mesh;
            s.young          = 1;
            s.parameterNames = {'m1','m2','a'};

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

        function createVolumeConstraint(obj)

            s.mesh         = obj.mesh;
            s.volumeTarget = 0.4;
            s.sideFields   = [1,2];

            obj.volume = LatticeVolumeConstraint(s);

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

            [mL,mU,aMax] = obj.createDesignParameters();

            limits = {[mL,mU],[mL,mU],[-aMax,aMax]};

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
            s.maxIter         = 1000;
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