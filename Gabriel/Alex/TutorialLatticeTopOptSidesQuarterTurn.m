classdef TutorialLatticeTopOptSidesQuarterTurn < handle
    properties (Access = public)
        complianceHistory
    end
    properties (Access = private)
        caseName
        networkFile
        young
        domainSize
        divisions
        clampLength
        loadLength
        totalLoad
        volumeTarget
        sideMin
        sideMax
        aMin
        aMax
        initialDensity
        initialA
        filterRadius
        etaStar
        delta
        tauMax
        maxIter
    end
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
        function obj = TutorialLatticeTopOptSidesQuarterTurn(changes)
            obj.init();
            if nargin > 0
                obj.applyChanges(changes);
            end
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
        end
        function compute(obj)
            tStart = tic;
            obj.createOptimizer();
            obj.complianceHistory = obj.compliance0*obj.cost.getIterationHistory();
            obj.saveMonitoring();
            obj.saveDesign();
            J = obj.complianceHistory;
            fprintf('%s: J(1) = %.4g, max J = %.4g, J(end) = %.4g, %d iterations, %.1f min\n',obj.caseName,J(1),max(J),J(end),numel(J)-1,toc(tStart)/60);
            fprintf('Load times displacement = %.4f\n',2*J(end));
        end
        function plotFields(obj)
            file = obj.getDesignFile();
            plotDesignFile(file,obj.caseName);
            plotDesignRhoBW(file,[obj.caseName,'_rho_bw.png']);
            plotWallDirection(file,QuarterTurnTransformation(),[],[obj.caseName,'_wallDirection']);
        end
        function plotCompliance(obj)
            J      = obj.loadComplianceHistory();
            result = struct('label',strrep(obj.caseName,'_',' '),'compliance',J(:),'final',J(end));
            plotComplianceHistories(result,[obj.caseName,'_compliance.png'],false,obj.volumeTarget);
        end
        function dehomogenizer = dehomogenize(obj,method,epsilon,elementsPerCell)
            if nargin < 4
                elementsPerCell = 16;
            end
            s.designFile      = obj.getDesignFile();
            s.epsilon         = epsilon;
            s.elementsPerCell = elementsPerCell;
            s.transformation  = QuarterTurnTransformation();
            s.parameterNames  = {'a'};
            switch method
                case 'simple'
                    dehomogenizer = LatticeDehomogenizer(s);
                case 'mapped'
                    dehomogenizer = MappedLatticeDehomogenizer(s);
                otherwise
                    error('TutorialLatticeTopOptSidesQuarterTurn:UnknownMethod','Unknown method: %s',method);
            end
            fileName = sprintf('%s_dehomo_%s_eps%s.png',obj.caseName,method,strrep(num2str(epsilon),'.','p'));
            plotLatticeBW(dehomogenizer,false,fileName);
            fprintf('%s, epsilon = %g: volume fraction %.4f, design %.4f, representable %.4f\n',method,epsilon, ...
                dehomogenizer.computeVolumeFraction(),dehomogenizer.computeTargetVolumeFraction(),dehomogenizer.computeRepresentableVolumeFraction());
        end
        function file = getDesignFile(obj)
            file = ['design_',obj.caseName,'.mat'];
        end
    end
    methods (Access = private)
        function init(obj)
            close all;
            obj.caseName       = 'perle_sidesQuarterTurn';
            obj.networkFile    = 'latticeSidesQuarterTurnNN_Adam_64x3_cholesky';
            obj.young          = 15;
            obj.domainSize     = [20,20];
            obj.divisions      = [80,80];
            obj.clampLength    = 10;
            obj.loadLength     = 1;
            obj.totalLoad      = -1;
            obj.volumeTarget   = 0.3;
            obj.sideMin        = 0.1;
            obj.sideMax        = 0.99;
            obj.aMin           = -2;
            obj.aMax           = 3;
            obj.initialDensity = 0.99;
            obj.initialA       = [0.5,0.5];
            obj.filterRadius   = [0,0,3];
            obj.etaStar        = 0.2;
            obj.delta          = 0.01;
            obj.tauMax         = 1000;
            obj.maxIter        = 1000;
        end
        function createMesh(obj)
            obj.mesh = TriangleMesh(obj.domainSize(1),obj.domainSize(2),obj.divisions(1),obj.divisions(2));
        end
        function createDesignVariable(obj)
            m0          = min(obj.sideMax,max(obj.sideMin,sqrt(1-obj.initialDensity)));
            sD.fields   = {obj.createUniformField(m0),obj.createUniformField(m0),obj.createInitialA()};
            sD.type     = 'MultiField';
            sD.plotting = true;
            obj.designVariable = DesignVariable.create(sD);
        end
        function f = createUniformField(obj,value)
            s.fHandle = @(x) value*ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            f         = AnalyticalFunction(s).project('P1');
        end
        function f = createInitialA(obj)
            aBottom   = obj.initialA(1);
            aTop      = obj.initialA(end);
            height    = obj.domainSize(2);
            s.fHandle = @(x) aBottom + (aTop-aBottom)*x(2,:,:)/height;
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            f         = AnalyticalFunction(s).project('P1');
        end
        function applyChanges(obj,changes)
            names = fieldnames(changes);
            for i = 1:numel(names)
                obj.(names{i}) = changes.(names{i});
            end
        end
        function createFilters(obj)
            h           = obj.mesh.computeMeanCellSize();
            obj.filters = cell(1,numel(obj.filterRadius));
            for iField = 1:numel(obj.filterRadius)
                obj.filters{iField} = obj.createFilter(obj.filterRadius(iField)*h);
            end
        end
        function f = createFilter(obj,epsilon)
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            if epsilon == 0
                s.filterType = 'LUMP';
                f = Filter.create(s);
            else
                s.filterType   = 'PDE';
                s.boundaryType = 'Neumann';
                s.metric       = 'Isotropy';
                f = Filter.create(s);
                f.updateEpsilon(epsilon);
            end
        end
        function createMaterial(obj)
            s.fileName       = obj.networkFile;
            s.mesh           = obj.mesh;
            s.young          = obj.young;
            s.parameterNames = {'m1','m2','a'};
            hm     = LatticeHomogenizedMaterialsReader(s);
            obj.C  = @(designVariables) hm.obtainTensor(designVariables);
            obj.dC = @(designVariables) hm.obtainTensorDerivatives(designVariables);
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
            obj.physicalProblem  = ElasticProblem(s);
        end
        function c = createComplianceFromConstitutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstitutiveTensor(s);
        end
        function s = createComplianceSettings(obj)
            s.mesh                       = obj.mesh;
            s.filters                    = obj.filters;
            s.complainceFromConstitutive = obj.createComplianceFromConstitutive();
            s.C                          = obj.C;
            s.dC                         = obj.dC;
        end
        function createCompliance(obj)
            obj.compliance = ComplianceFunctional(obj.createComplianceSettings());
        end
        function computeInitialCompliance(obj)
            s               = obj.createComplianceSettings();
            s.value0        = 1;
            physical        = ComplianceFunctional(s);
            obj.compliance0 = physical.computeFunctionAndGradient(obj.designVariable);
        end
        function createVolumeConstraint(obj)
            s.mesh         = obj.mesh;
            s.volumeTarget = obj.volumeTarget;
            s.sideFields   = [1,2];
            obj.volume     = LatticeVolumeConstraint(s);
        end
        function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = RecordingCost(Cost(s));
        end
        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end
        function M = createMassMatrix(obj)
            test    = LagrangianFunction.create(obj.mesh,1,'P1');
            trial   = LagrangianFunction.create(obj.mesh,1,'P1');
            MSingle = IntegrateLHS(@(u,v) DP(v,u),test,trial,obj.mesh,'Domain');
            MSingle = diag(sum(MSingle,1));
            nFields = obj.designVariable.obtainNumberVariables();
            M       = kron(speye(nFields),MSingle);
        end
        function limits = createFieldLimits(obj)
            limits = {[obj.sideMin,obj.sideMax],[obj.sideMin,obj.sideMax],[obj.aMin,obj.aMax]};
        end
        function createPrimalUpdater(obj)
            fields = obj.designVariable.obtainDomainFunction();
            limits = obj.createFieldLimits();
            s.lb   = [];
            s.ub   = [];
            for iField = 1:numel(fields)
                n    = numel(fields{iField}.fValues);
                s.lb = [s.lb; limits{iField}(1)*ones(n,1)];
                s.ub = [s.ub; limits{iField}(2)*ones(n,1)];
            end
            s.tauMax = obj.tauMax;
            s.tau    = [];
            obj.primalUpdater = ProjectedGradient(s);
        end
        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = obj.maxIter;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.primal         = 'PROJECTED GRADIENT';
            s.delta          = obj.delta;
            s.etaStar        = obj.etaStar;
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
            xMax      = max(obj.mesh.coord(:,1));
            yMid      = max(obj.mesh.coord(:,2))/2;
            tol       = 1e-9;
            yBoundary = unique(obj.mesh.coord(obj.mesh.coord(:,1) == xMax,2));
            spacing   = min(diff(yBoundary));
            traction  = obj.totalLoad/obj.loadLength;
            isClamp        = @(coor) coor(:,1) == 0 & abs(coor(:,2)-yMid) <= obj.clampLength/2 + tol;
            isLoadInterior = @(coor) coor(:,1) == xMax & abs(coor(:,2)-yMid) < obj.loadLength/2 - tol;
            isLoadEnd      = @(coor) coor(:,1) == xMax & abs(abs(coor(:,2)-yMid) - obj.loadLength/2) <= tol;
            sDir.domain    = isClamp;
            sDir.direction = [1,2];
            sDir.value     = 0;
            sPL{1}.domain    = isLoadInterior;
            sPL{1}.direction = 2;
            sPL{1}.value     = traction*spacing;
            sPL{2}.domain    = isLoadEnd;
            sPL{2}.direction = 2;
            sPL{2}.value     = traction*spacing/2;
            s.dirichletFun = DirichletCondition(obj.mesh,sDir);
            s.pointloadFun = [];
            for i = 1:numel(sPL)
                s.pointloadFun = [s.pointloadFun,TractionLoad(obj.mesh,sPL{i},'DIRAC')];
            end
            s.periodicFun = [];
            s.mesh        = obj.mesh;
            bc = BoundaryConditions(s);
        end
        function saveDesign(obj)
            fields        = obj.designVariable.fields;
            design.coord  = obj.mesh.coord;
            design.connec = obj.mesh.connec;
            design.m1     = obj.filters{1}.compute(fields{1},2).fValues;
            design.m2     = obj.filters{2}.compute(fields{2},2).fValues;
            design.a      = obj.filters{3}.compute(fields{3},2).fValues;
            design.rho    = 1 - design.m1.*design.m2;
            design.xi     = atan2(design.m1,design.m2);
            design.mU     = obj.sideMax;
            design.compliance = obj.complianceHistory;
            save(obj.getDesignFile(),'-struct','design');
        end
        function saveMonitoring(obj)
            figs = findall(0,'Type','figure');
            for k = 1:numel(figs)
                print(figs(k),sprintf('%s_monitoring_fig%d.png',obj.caseName,figs(k).Number),'-dpng','-r200');
            end
        end
        function J = loadComplianceHistory(obj)
            if isempty(obj.complianceHistory)
                d = load(obj.getDesignFile(),'compliance');
                obj.complianceHistory = d.compliance;
            end
            J = obj.complianceHistory;
        end
    end
end