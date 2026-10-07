classdef TutorialFirstTwovariable < handle

    properties (Access = private)
        mesh
        filterPDE
        filterPDERho
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
        % penalty
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
        function hist = getOptimizationDiagnostics(obj)

            hist.cost       = obj.optimizer.costHistory;
            hist.volume     = obj.optimizer.volumeConstraintHistory;

            hist.normDJb    = obj.optimizer.normDJbHistory;
            hist.normDJrho  = obj.optimizer.normDJrhoHistory;

            hist.normDxJb   = obj.optimizer.normDxJbHistory;
            hist.normDxJrho = obj.optimizer.normDxJrhoHistory;

            hist.normMeritB   = obj.optimizer.normMeritBHistory;
            hist.normMeritRho = obj.optimizer.normMeritRhoHistory;

            hist.deltaB   = obj.optimizer.deltaBHistory;
            hist.deltaRho = obj.optimizer.deltaRhoHistory;

            hist.fracBLower = obj.optimizer.fracBLowerHistory;
            hist.fracBUpper = obj.optimizer.fracBUpperHistory;

            hist.fracRhoLower = obj.optimizer.fracRhoLowerHistory;
            hist.fracRhoUpper = obj.optimizer.fracRhoUpperHistory;
            hist.fracRhoBound = obj.optimizer.fracRhoBoundHistory;

            hist.meanB   = obj.optimizer.meanBHistory;
            hist.meanRho = obj.optimizer.meanRhoHistory;

            hist.minB = obj.optimizer.minBHistory;
            hist.maxB = obj.optimizer.maxBHistory;

            hist.minRho = obj.optimizer.minRhoHistory;
            hist.maxRho = obj.optimizer.maxRhoHistory;
            hist.tau = obj.optimizer.tauHistory;

            hist.eta = obj.optimizer.etaHistory;

            hist.etaMax = obj.optimizer.etaMaxHistory;

            hist.lineSearchTrials = ...
                obj.optimizer.lineSearchTrialsHistory;

        end
        function data = getFinalKKTDiagnostics(obj)
            data = obj.optimizer.getFinalKKTDiagnostics();
        end

        function hist = getComplianceHistory(obj)
            hist = obj.optimizer.complianceHistory;
        end

        function hist = getVolumeConstraintHistory(obj)
            hist = obj.optimizer.volumeConstraintHistory;
        end
        function b = getBValues(obj)
            
            b = obj.designVariable.funB.fValues;
        end

        function rho = getRhoValues(obj)
            
            rho = obj.designVariable.funRho.fValues;
        end

        function bSmooth = getBSmoothValues(obj)
           
            bSmooth = obj.bSmooth.fValues;
        end

        function rhoSmooth = getRhoSmoothValues(obj)
           
            rhoSmooth = obj.rhoSmooth.fValues;
        end

        function fields = getDesignVariableValues(obj)
            

            fields.b   = obj.designVariable.funB.fValues;
            fields.rho = obj.designVariable.funRho.fValues;

            fields.bSmooth   = obj.bSmooth.fValues;
            fields.rhoSmooth = obj.rhoSmooth.fValues;
        end

        function meshData = getMeshData(obj)
            

            meshData.coord  = obj.mesh.coord;
            meshData.connec = obj.mesh.connec;
            meshData.ndim   = obj.mesh.ndim;
            meshData.nnodes = obj.mesh.nnodes;
        end

        function obj = TutorialFirstTwovariable()
            obj.init();
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilters();
            obj.createBaseDomain();
            % obj.createInterfacePenalty();
            obj.createMaterial();
            obj.createElasticProblem();
            obj.createComplianceFromConstitutive();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createOptimizer();
            obj.bSmooth   = obj.filterPDE.compute(obj.designVariable.funB,   3);
            obj.rhoSmooth = obj.filterPDERho.compute(obj.designVariable.funRho, 3);
           
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            obj.mesh = TriangleMesh(2, 1, 150, 75);
        end

        function createDesignVariable(obj)
            s_b.fHandle = @(x) zeros(size(x(1,:,:)));
            s_b.ndimf   = 1;
            s_b.mesh    = obj.mesh;
            aFunB = AnalyticalFunction(s_b);
            funB  = aFunB.project('P1');

            s_rho.fHandle = @(x) 0.998*ones(size(x(1,:,:)));
            s_rho.ndimf   = 1;
            s_rho.mesh    = obj.mesh;
            aFunRho = AnalyticalFunction(s_rho);
            funRho  = aFunRho.project('P1');

            sD.funB     = funB;
            sD.fun      = funRho;
            sD.mesh     = obj.mesh;
            sD.type     = 'MicroWithDensity';
            sD.plotting = true;
            obj.designVariable = DesignVariable.create(sD);
        end

        function createFilters(obj)

            h = obj.mesh.computeMeanCellSize();

            % ========================================================
            % FILTER FOR b
            % ========================================================

            sB.filterType    = 'PDE';
            sB.mesh          = obj.mesh;
            sB.boundaryType  = 'Neumann';
            sB.metric        = 'Isotropy';
            sB.trial         = ...
                LagrangianFunction.create(obj.mesh,1,'P1');

            fB = Filter.create(sB);

            epsilonB = 6*h;

            fB.updateEpsilon(epsilonB);

            obj.filterPDE = fB;


            % ========================================================
            % FILTER FOR rho
            % ========================================================

            sR.filterType    = 'PDE';
            sR.mesh          = obj.mesh;
            sR.boundaryType  = 'Neumann';
            sR.metric        = 'Isotropy';
            sR.trial         = ...
                LagrangianFunction.create(obj.mesh,1,'P1');

            fR = Filter.create(sR);

            epsilonRho = 3*h;

            fR.updateEpsilon(epsilonRho);

            obj.filterPDERho = fR;


            % ========================================================
            % LUMP FILTER - keep available
            % ========================================================

            sL.filterType = 'LUMP';
            sL.mesh       = obj.mesh;
            sL.trial      = ...
                LagrangianFunction.create(obj.mesh,1,'P1');

            obj.filterLUMP = Filter.create(sL);

        end

        % function createInterfacePenalty(obj)
        %     sF.mesh       = obj.mesh;
        %     sF.filterType = 'PDE';
        %     sF.trial      = ...
        %         LagrangianFunction.create( ...
        %         obj.mesh,1,'P1');
        % 
        %     f = Filter.create(sF);
        % 
        %     h = obj.mesh.computeMeanCellSize();
        % 
        %     s.mesh        = obj.mesh;
        %     s.uMesh       = obj.baseDomain;
        %     s.filter      = f;
        % 
        %     s.epsilon     = 3*h;
        % 
        %     s.value0      = 8;
        % 
        %     s.signInitial = -3;
        %     s.signFinal   = -1.2;
        % 
        %     s.tarVolume   = 0.4;
        % 
        %     obj.penalty = ...
        %         InterfaceFunctionalMicroDensity(s);
        % 
        % end

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
            s.fileName    = 'Homogenizationtwovariables24NewT';
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
            s.filterDensity              = obj.filterPDERho;
            s.complainceFromConstitutive = obj.compliance;
            s.material                   = obj.materialMicro;
            c = ComplianceFunctionalMicroDensity(s);
        end

        function createVolumeConstraint(obj)
            s.mesh         = obj.mesh;
            s.filter       = obj.filterPDERho;
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

            % s.shapeFunctions{1} = ...
            %     obj.createCompliance();
            % 
            % s.shapeFunctions{2} = ...
            %     obj.penalty;
            % 
            % s.weights = [1,1];
            % 
            % s.Msmooth = ...
            %     obj.createMassMatrix();
            % 
            % obj.cost = Cost(s);

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
            MSingle = diag(sum(MSingle,1));
            n = test.nDofs;
            Z = sparse(n, n);
            M = [MSingle, Z; Z, MSingle];
        end

        function p = createPrimalUpdater(obj)
            n    = obj.mesh.nnodes;
            s.lb = [-ones(n,1);  1e-6*ones(n,1)];
            s.ub = [ ones(n,1);  0.998*ones(n,1)];
            s.tauMax = 100;
            s.tau    = [];
            p = ProjectedGradient(s);
        end

        function createOptimizer(obj)
            s.monitoring      = true;
            s.cost            = obj.cost;
            s.constraint      = obj.constraint;
            s.designVariable  = obj.designVariable;
            s.maxIter         = 1500;
            s.tolerance       = 1e-8;
            s.constraintCase  = {'EQUALITY'};
            s.etaNorm         = 0.01;
            s.gJFlowRatio     = 0.8;
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