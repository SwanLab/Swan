classdef TopOptPunzonLevelSet < handle

    properties (Access = private)
        filename
        mesh
        filter
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        volume
        cost
        constraint
        primalUpdater
        % dualVariable
        optimizer
    end

    methods (Access = public)

        function obj = TopOptPunzonLevelSet()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            % obj.createDualVariable();
            obj.createOptimizer();

            figure(2)
            saveas(gcf,'Monitoring_PunzonLevelSet.fig');
            obj.designVariable.fun.print('fValues_PunzonLevelSet');
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            file = 'punzon';
            obj.filename = file;
            a.fileName = file;
            s = FemDataContainer(a);
            obj.mesh = s.mesh;
        end

        % function createDesignVariable(obj)
        %     s.fHandle = @(x) ones(size(x(1,:,:)));
        %     s.ndimf   = 1;
        %     s.mesh    = obj.mesh;
        %     aFun      = AnalyticalFunction(s);
        %     s.fun     = aFun.project('P1');
        %     s.fun.setFValues(0.4*ones(size(s.fun.fValues,1),1));
        %     s.mesh    = obj.mesh;
        %     s.type = 'Density';
        %     s.plotting = true;
        %     dens    = DesignVariable.create(s);
        %     obj.designVariable = dens;
        % end

        function createDesignVariable(obj)
            s.type = 'Full';
            g      = GeometricalFunction(s);
            lsFun  = g.computeLevelSetFunction(obj.mesh);
            s.fun  = lsFun;
            s.mesh = obj.mesh;
            s.type = 'LevelSet';
            s.plotting = true;
            ls     = DesignVariable.create(s);
            obj.designVariable = ls;
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
            s.dim            = '3D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;
        end

        function m = createMaterial(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = f{1}.project('P1');            
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '3D';
            s.mesh                 = obj.mesh;
            m = Material.create(s);
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '3D';
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
            s.material                    = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function uMesh = createBaseDomain(obj)
            sG.type          = 'Full';
            sG.length        = 1;
            sG.xCoorCenter   = 1.5;
            sG.yCoorCenter   = 0.5;
            g                = GeometricalFunction(sG);
            lsFun            = g.computeLevelSetFunction(obj.mesh);
            levelSet         = lsFun.fValues;
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh            = UnfittedMesh(s);
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

        % function createDualVariable(obj)
        %     s.nConstraints   = 1;
        %     l                = DualVariable(s);
        %     obj.dualVariable = l;
        % end

        function createPrimalUpdater(obj)
            s.mesh = obj.mesh;
            obj.primalUpdater = SLERP(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 15;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.primalUpdater  = obj.primalUpdater;
            s.etaNorm        = 0.1;
            s.etaNormMin     = 0.05;
            s.gJFlowRatio    = 0.05;
            s.etaMax         = 1;
            s.etaMaxMin      = 0.1;
            %s.type           = '0';
            s.gif = false;
            s.gifName = [];
            s.printing = true;
            s.printName = 'PunzonLevelSet';
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj)
            zMin = min(obj.mesh.coord(:,3));
            zMax = max(obj.mesh.coord(:,3));

            isTop = @(coor) abs(coor(:,3) - zMax) < 1e-6; % tornillos
            isForce = @(x) (x(3,:,:) - zMin) < 1e-2;

            sDir{1}.domain    = @(coor) isTop(coor);
            sDir{1}.direction = [1,2,3];
            sDir{1}.value     = 0;

            sPL{1}.domain    = isForce;
            [bMesh,~]=obj.mesh.createSingleBoundaryMesh();
            sPL{1}.fun       = ConstantFunction.create([0,0,1],bMesh);


            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = TractionLoad(obj.mesh, sPL{i}, 'FUNCTION');
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun = [];
            s.mesh        = obj.mesh;

            bc = BoundaryConditions(s);
        end
    end

    methods (Static, Access=private)
        function sCond = computeCondition(conditions)
            nodes = @(coor) 1:size(coor,1);
            dirs  = unique(conditions(:,2));
            j     = 0;
            for k = 1:length(dirs)
                rowsDirk = ismember(conditions(:,2),dirs(k));
                u        = unique(conditions(rowsDirk,3));
                for i = 1:length(u)
                    rows   = conditions(:,3)==u(i) & rowsDirk;
                    isCond = @(coor) ismember(nodes(coor),conditions(rows,1));
                    j      = j+1;
                    sCond{j}.domain    = @(coor) isCond(coor);
                    sCond{j}.direction = dirs(k);
                    sCond{j}.value     = u(i);
                end
            end
        end
    end
end