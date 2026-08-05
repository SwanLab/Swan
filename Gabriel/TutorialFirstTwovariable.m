classdef TutorialFirstTwovariable < handle

    properties (Access = private)
        mesh
        filterPDE
        filterLUMP
        % designVariable
        physicalProblem
        materialMicro
        compliance
        volume
        cost
        constraint
        optimizer
        baseDomain
        initialB
        initialBCoord
    end
    properties (Access = public)
        designVariable
        bSmooth
        rhoSmooth
    end

    methods (Access = public)
        function hist = getCostHistory(obj)
            hist = obj.optimizer.costHistory;
        end

        function hist = getComplianceHistory(obj)
            hist = obj.optimizer.complianceHistory;
        end

        function hist = getVolumeConstraintHistory(obj)
            hist = obj.optimizer.volumeConstraintHistory;
        end
        function b = getBValues(obj)
            % Campo b final antes do filtro
            b = obj.designVariable.funB.fValues;
        end

        function rho = getRhoValues(obj)
            % Campo rho final antes do filtro
            rho = obj.designVariable.funRho.fValues;
        end

        function bSmooth = getBSmoothValues(obj)
            % Campo b final depois do filtro PDE
            bSmooth = obj.bSmooth.fValues;
        end

        function rhoSmooth = getRhoSmoothValues(obj)
            % Campo rho final depois do filtro LUMP
            rhoSmooth = obj.rhoSmooth.fValues;
        end

        function fields = getDesignVariableValues(obj)
            % Estrutura unica para o script comparativo

            fields.b   = obj.designVariable.funB.fValues;
            fields.rho = obj.designVariable.funRho.fValues;

            fields.bSmooth   = obj.bSmooth.fValues;
            fields.rhoSmooth = obj.rhoSmooth.fValues;
        end

        function meshData = getMeshData(obj)
            % Dados numericos da malha

            meshData.coord  = obj.mesh.coord;
            meshData.connec = obj.mesh.connec;
            meshData.ndim   = obj.mesh.ndim;
            meshData.nnodes = obj.mesh.nnodes;
        end

        function obj = TutorialFirstTwovariable(bInitial,coordInitial)

            if nargin < 1
                bInitial = [];
            end

            if nargin < 2
                coordInitial = [];
            end

            obj.initialB = bInitial;
            obj.initialBCoord = coordInitial;

            obj.init();
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilters();
            obj.createBaseDomain();
            obj.createMaterial();
            obj.createElasticProblem();
            obj.createComplianceFromConstitutive();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createOptimizer();

            obj.bSmooth = obj.filterPDE.compute( ...
                obj.designVariable.funB,3);

            obj.rhoSmooth = obj.filterLUMP.compute( ...
                obj.designVariable.funRho,2);
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            obj.mesh = TriangleMesh(2, 1, 120, 100);
        end

        function createDesignVariable(obj)

            %--------------------------------------------------------------
            % Campo inicial de b
            %--------------------------------------------------------------
            if isempty(obj.initialB)

                s_b.fHandle = @(x) zeros(size(x(1,:,:)));

                fprintf('\nInitial guess de b uniforme: b0 = 0.\n');

            else

                bInitial     = obj.initialB(:);
                coordInitial = obj.initialBCoord;

                if isempty(coordInitial)
                    error(['As coordenadas da malha do campo inicial de b ', ...
                        'tambem precisam ser fornecidas.']);
                end

                if size(coordInitial,1) ~= numel(bInitial)
                    error(['O numero de coordenadas (%d) nao coincide com ', ...
                        'o numero de valores de b (%d).'], ...
                        size(coordInitial,1),numel(bInitial));
                end

                if size(coordInitial,2) ~= 2
                    error('coordInitial deve possuir duas colunas: x e y.');
                end

                if any(~isfinite(bInitial))
                    error('O campo inicial de b possui NaN ou Inf.');
                end

                if any(bInitial < -0.8) || any(bInitial > 0.8)
                    error(['O campo inicial de b esta fora dos limites ', ...
                        '[-0.8,0.8]. Min = %.6f, Max = %.6f.'], ...
                        min(bInitial),max(bInitial));
                end

                % Verifica se a malha é exatamente a mesma
                sameSize = size(coordInitial,1) == obj.mesh.nnodes;

                if sameSize
                    maxCoordDifference = max(abs( ...
                        coordInitial(:)-obj.mesh.coord(:)));
                else
                    maxCoordDifference = Inf;
                end

                fprintf('\nInitial guess de b carregado.\n');
                fprintf('Numero de valores = %d\n',numel(bInitial));
                fprintf('b inicial min/max = %.6f  %.6f\n', ...
                    min(bInitial),max(bInitial));
                fprintf('Diferenca maxima entre malhas = %.6e\n', ...
                    maxCoordDifference);

                % Interpolador espacial do campo salvo
                Fb = scatteredInterpolant( ...
                    coordInitial(:,1), ...
                    coordInitial(:,2), ...
                    bInitial, ...
                    'linear', ...
                    'nearest');

                s_b.fHandle = @(x) Fb(x(1,:,:),x(2,:,:));
            end

            s_b.ndimf = 1;
            s_b.mesh  = obj.mesh;

            aFunB = AnalyticalFunction(s_b);
            funB  = aFunB.project('P1');

            %--------------------------------------------------------------
            % Campo inicial de rho
            %--------------------------------------------------------------
            s_rho.fHandle = @(x) 0.97*ones(size(x(1,:,:)));
            s_rho.ndimf   = 1;
            s_rho.mesh    = obj.mesh;

            aFunRho = AnalyticalFunction(s_rho);
            funRho  = aFunRho.project('P1');

            %--------------------------------------------------------------
            % Variável conjunta [b; rho]
            %--------------------------------------------------------------
            sD.funB     = funB;
            sD.fun      = funRho;
            sD.mesh     = obj.mesh;
            sD.type     = 'MicroWithDensity';
            sD.plotting = true;

            obj.designVariable = DesignVariable.create(sD);
        end

        function createFilters(obj)
            eOverhmin = 6;
            s.filterType   = 'PDE';
            s.mesh         = obj.mesh;
            s.boundaryType = 'Neumann';
            s.metric       = 'Isotropy';
            s.trial        = LagrangianFunction.create(obj.mesh, 1, 'P1');
            f = Filter.create(s);
            epsilon = eOverhmin * obj.mesh.computeMeanCellSize();
            f.updateEpsilon(epsilon);
            obj.filterPDE = f;

            sL.filterType = 'LUMP';
            sL.mesh       = obj.mesh;
            sL.trial      = LagrangianFunction.create(obj.mesh, 1, 'P1');
            obj.filterLUMP = Filter.create(sL);
        end

        function createBaseDomain(obj)
            levelSet = -ones(obj.mesh.nnodes, 1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
            uMesh.compute(levelSet);
            obj.baseDomain = uMesh;
        end

        function createMaterial(obj)
            s.type        = 'HomogenizedMicroDensityFixed';
            s.mesh        = obj.mesh;
            s.young       = 1.0;
            s.fileName    = 'Homogenizationtwovariables15';
            s.density     = obj.designVariable;
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

        function createComplianceFromConstitutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            obj.compliance = ComplianceFromConstitutiveTensor(s);
        end

        function c = createCompliance(obj)
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filterPDE;
            s.filterDensity              = obj.filterLUMP;
            s.complainceFromConstitutive = obj.compliance;
            s.material                   = obj.materialMicro;
            c = ComplianceFunctionalMicroDensity(s);
        end

        function createVolumeConstraint(obj)
            s.mesh         = obj.mesh;
            s.filter       = obj.filterLUMP;
            s.test         = LagrangianFunction.create(obj.mesh, 1, 'P1');
            s.volumeTarget = 0.4;
            s.uMesh        = obj.baseDomain;
            obj.volume     = VolumeConstraintMicroDensity(s);
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.createCompliance();
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function M = createMassMatrix(obj)
            test  = LagrangianFunction.create(obj.mesh, 1, 'P1');
            trial = LagrangianFunction.create(obj.mesh, 1, 'P1');
            MSingle = IntegrateLHS(@(u,v) DP(v,u), test, trial, obj.mesh, 'Domain');
            n = test.nDofs;
            Z = sparse(n, n);
            M = [MSingle, Z; Z, MSingle];
        end

        function p = createPrimalUpdater(obj)
            n    = obj.mesh.nnodes;
            s.lb = [-0.8*ones(n,1);  1e-6*ones(n,1)];
            s.ub = [ 0.8*ones(n,1);  0.97*ones(n,1)];
            s.tauMax = 500;
            s.tau    = [];
            p = ProjectedGradient(s);
        end

        function createOptimizer(obj)
            s.monitoring      = true;
            s.cost            = obj.cost;
            s.constraint      = obj.constraint;
            s.designVariable  = obj.designVariable;
            s.maxIter         = 3000;
            s.tolerance       = 1e-8;
            s.constraintCase  = {'EQUALITY'};
            s.etaNorm         = 0.01;
            s.gJFlowRatio     = 0.2;
            s.primalUpdater   = obj.createPrimalUpdater();
            s.gif             = false;
            s.gifName         = [];
            s.printing        = false;
            s.printName       = [];
            
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
            
            
        end

        function bc = createBoundaryConditions(obj)
            xMin = min(obj.mesh.coord(:,1));
            xMax = max(obj.mesh.coord(:,1));
            yMax = max(obj.mesh.coord(:,2));

            isDir   = @(coor) abs(coor(:,1) - xMin) < 1e-12;
            isForce = @(coor) abs(coor(:,1) - xMax) < 1e-12 & ...
                              coor(:,2) >= 0.4*yMax & ...
                              coor(:,2) <= 0.6*yMax;

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1, 2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = -1;

            dirichletFun = [];
            for i = 1:numel(sDir)
                dirichletFun = [dirichletFun, DirichletCondition(obj.mesh, sDir{i})];
            end

            pointloadFun = [];
            for i = 1:numel(sPL)
                pointloadFun = [pointloadFun, TractionLoad(obj.mesh, sPL{i}, 'DIRAC')];
            end

            s.dirichletFun = dirichletFun;
            s.pointloadFun = pointloadFun;
            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
        end

    end

end