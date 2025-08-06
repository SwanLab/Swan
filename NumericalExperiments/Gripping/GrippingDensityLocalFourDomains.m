classdef GrippingDensityLocalFourDomains < handle

    properties (Access = private)
        mesh
        filename
        filterCompliance
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        volume
        perimeterLU
        perimeterLD
        perimeterRU
        perimeterRD
        cost
        constraint
        primalUpdater
        optimizer
    end

    methods (Access = public)

        function obj = GrippingDensityLocalFourDomains()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilterCompliance();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createNonSelfAdjCompliance();
            obj.createVolumeConstraint();
            obj.createPerimeterLU();
            obj.createPerimeterLD();
            obj.createPerimeterRU();
            obj.createPerimeterRD();
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();

            d = obj.designVariable;
            saveas(gcf,'NumericalExperiments/Gripping/MonitDensityFourDomains.fig');
            save('NumericalExperiments/Gripping/DesVarDensityFourDomains.mat','d');
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            file = 'Gripping';
            obj.filename = file;
            a.fileName = file;
            s = FemDataContainer(a);
            obj.mesh = s.mesh;        
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);
            
            sD.fun      = aFun.project('P1');
            sD.mesh     = obj.mesh;
            sD.type     = 'Density';
            sD.plotting = false;
            dens        = DesignVariable.create(sD);
            obj.designVariable = dens;
        end

        function createFilterCompliance(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filterCompliance = f;
        end

        function f = createFilterPerimeter(obj)
            s.trial   = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh    = obj.mesh;
            s.LHStype = 'StiffnessMass';
            f         = FilterPDE(s);
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
            f = obj.designVariable.fun;           
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            s.mesh                 = obj.mesh;
            m = Material.create(s);
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function createNonSelfAdjCompliance(obj)
            s.mesh         = obj.mesh;
            s.filter       = obj.filterCompliance;
            s.material     = obj.createMaterial();
            s.stateProblem = obj.physicalProblem;
            s.filename     = obj.filename;
            c = NonSelfAdjointComplianceFunctional(s);
            obj.compliance = c;
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filterCompliance;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.6;
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function base = createBaseDomain(obj,x0,y0)
            s.type             = 'Rectangle';
            s.xCoorCenter      = x0;
            s.yCoorCenter      = y0;
            s.xSide            = 0.5;
            s.ySide            = 0.5;
            g                  = GeometricalFunction(s);
            lsFun              = g.computeLevelSetFunction(obj.mesh);
            sUm.backgroundMesh = obj.mesh;
            sUm.boundaryMesh   = obj.mesh.createBoundaryMesh();
            base               = UnfittedMesh(sUm);
            base.compute(lsFun.fValues);
        end

        function createPerimeterLU(obj)
            s.mesh         = obj.mesh;
            s.uMesh        = obj.createBaseDomain(0.25,0.75);
            s.epsilon      = 5*obj.mesh.computeMeanCellSize();
            s.filter       = obj.createFilterPerimeter();
            s.value0       = 1;
            s.minEpsilon   = 1.5*obj.mesh.computeMeanCellSize();
            s.target       = 1.35*(3/4);
            obj.perimeterLU = PerimeterConstraint(s);
        end

        function createPerimeterLD(obj)
            s.mesh         = obj.mesh;
            s.uMesh        = obj.createBaseDomain(0.25,0.25);
            s.epsilon      = 5*obj.mesh.computeMeanCellSize();
            s.filter       = obj.createFilterPerimeter();
            s.value0       = 1;
            s.minEpsilon   = 1.5*obj.mesh.computeMeanCellSize();
            s.target       = 1.35*(3/4);
            obj.perimeterLD = PerimeterConstraint(s);
        end

        function createPerimeterRU(obj)
            s.mesh         = obj.mesh;
            s.uMesh        = obj.createBaseDomain(0.75,0.75);
            s.epsilon      = 5*obj.mesh.computeMeanCellSize();
            s.filter       = obj.createFilterPerimeter();
            s.value0       = 1;
            s.minEpsilon   = 1.5*obj.mesh.computeMeanCellSize();
            s.target       = 1.06*(3/4);
            obj.perimeterRU = PerimeterConstraint(s);
        end

        function createPerimeterRD(obj)
            s.mesh         = obj.mesh;
            s.uMesh        = obj.createBaseDomain(0.75,0.25);
            s.epsilon      = 5*obj.mesh.computeMeanCellSize();
            s.filter       = obj.createFilterPerimeter();
            s.value0       = 1;
            s.minEpsilon   = 1.5*obj.mesh.computeMeanCellSize();
            s.target       = 1.06*(3/4);
            obj.perimeterRD = PerimeterConstraint(s);
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            n = obj.mesh.nnodes;
            h = obj.mesh.computeMinCellSize();
            M = h^2*sparse(1:n,1:n,ones(1,n),n,n);
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.shapeFunctions{2} = obj.perimeterLU;
            s.shapeFunctions{3} = obj.perimeterLD;
            s.shapeFunctions{4} = obj.perimeterRU;
            s.shapeFunctions{5} = obj.perimeterRD;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createPrimalUpdater(obj)
            s.ub     = 1;
            s.lb     = 0;
            s.tauMax = 1000;
            obj.primalUpdater = ProjectedGradient(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 2000;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY','INEQUALITY','INEQUALITY','INEQUALITY','INEQUALITY'};
            s.primal         = 'PROJECTED GRADIENT';
            s.etaNorm        = 0.02;
            s.gJFlowRatio    = 0.2;
            s.primalUpdater  = obj.primalUpdater;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj)
            femReader = FemInputReaderGiD();
            s         = femReader.read(obj.filename);
            sPL       = obj.computeCondition(s.pointload);
            sDir      = obj.computeCondition(s.dirichlet);

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
