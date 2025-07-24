classdef GripperProblemLevelSetLocalPerimeter < handle

    properties (Access = private)
        filename
        filter
        materialInterpolator
        physicalProblem
        compliance
        volume
        globalPerimeter
        hinge1Perimeter
        hinge2Perimeter
        hinge3Perimeter
        hinge4Perimeter
        cost
        constraint
        dualVariable
        optimizer
    end

    properties (Access = public)
        mesh
        designVariable
    end

    methods (Access = public)

        function obj = GripperProblemLevelSetLocalPerimeter(r,gJ,w)
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();
            obj.createNonSelfAdjCompliance();
            obj.createVolumeConstraint();
            obj.createPerimeterConstraintHinge1(r);
            obj.createPerimeterConstraintHinge2(r);
            obj.createPerimeterConstraintHinge3(r);
            obj.createPerimeterConstraintHinge4(r);
            obj.createGlobalPerimeterConstraint();
            obj.createCost(w);
            obj.createConstraint();
            obj.createDualVariable();
            obj.createOptimizer(gJ);

            fileLocation = 'C:\Users\Biel\Desktop\UNI\TFG\ResultatsNormP_Density\00. From Batch';
            
            vtuName = fullfile(fileLocation, sprintf('Topology_Gripper_localPerimeter_r%.2f_gJ%.2f_eta0.02_w%.2f_LevelSet',r,gJ,w));
            obj.designVariable.fun.print(vtuName);
            
            figure(2)
            set(gcf, 'Position', get(0, 'Screensize'));
            fileName1 = fullfile(fileLocation, sprintf('Monitoring_Gripper_localPerimeter_r%.2f_gJ%.2f_eta0.02_w%.2f_LevelSet.fig',r,gJ,w));
            fileName2 = fullfile(fileLocation, sprintf('Monitoring_Gripper_localPerimeter_r%.2f_gJ%.2f_eta0.02_w%.2f_LevelSet.png',r,gJ,w));
            savefig(fileName1);
            print(fileName2,'-dpng','-r300');
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
            E0   = 1e-3;
            nu0  = 1/3;
            E1   = 1;
            nu1  = 1/3;
            ndim = 2;

            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.typeOfMaterial = 'ISOTROPIC';
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

        function c = createComplianceFromConstiutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstitutiveTensor(s);
        end

        function createNonSelfAdjCompliance(obj)
            s.mesh         = obj.mesh;
            s.filter       = obj.filter;
            s.material     = obj.createMaterial();
            s.stateProblem = obj.physicalProblem;
            s.filename     = obj.filename;
            c = NonSelfAdjointComplianceFunctional(s);
            obj.compliance = c;
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.6;
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createPerimeterConstraintHinge1(obj,r)
            s.mesh              = obj.mesh;
            s.uMesh             = obj.computeHinge1Domain(r);
            s.filter            = obj.createFilterPerimeter();
            s.epsilon           = 4*obj.mesh.computeMeanCellSize();
            s.minEpsilon        = 4*obj.mesh.computeMeanCellSize();
            s.value0            = 0.01;
            s.target            = 0.001;
            obj.hinge1Perimeter = LocalPerimeterConstraint(s);
        end

        function createPerimeterConstraintHinge2(obj,r)
            s.mesh              = obj.mesh;
            s.uMesh             = obj.computeHinge2Domain(r);
            s.filter            = obj.createFilterPerimeter();
            s.epsilon           = 4*obj.mesh.computeMeanCellSize();
            s.minEpsilon        = 4*obj.mesh.computeMeanCellSize();
            s.value0            = 0.01;
            s.target            = 0.001;
            obj.hinge2Perimeter = LocalPerimeterConstraint(s);
        end

        function createPerimeterConstraintHinge3(obj,r)
            s.mesh              = obj.mesh;
            s.uMesh             = obj.computeHinge3Domain(r);
            s.filter            = obj.createFilterPerimeter();
            s.epsilon           = 4*obj.mesh.computeMeanCellSize();
            s.minEpsilon        = 4*obj.mesh.computeMeanCellSize();
            s.value0            = 0.01;
            s.target            = 0.001;
            obj.hinge3Perimeter = LocalPerimeterConstraint(s);
        end

        function createPerimeterConstraintHinge4(obj,r)
            s.mesh              = obj.mesh;
            s.uMesh             = obj.computeHinge4Domain(r);
            s.filter            = obj.createFilterPerimeter();
            s.epsilon           = 4*obj.mesh.computeMeanCellSize();
            s.minEpsilon        = 4*obj.mesh.computeMeanCellSize();
            s.value0            = 0.01;
            s.target            = 0.001;
            obj.hinge4Perimeter = LocalPerimeterConstraint(s);
        end

        function uMesh = computeHinge1Domain(obj,r)
            s.type        = 'Circle';
            s.radius      = r*obj.mesh.computeMeanCellSize();
            s.xCoorCenter = 0.5;
            s.yCoorCenter = 0.75;
            g             = GeometricalFunction(s);
            lsFun         = g.computeLevelSetFunction(obj.mesh);
            levelSet      = lsFun.fValues;
            sUm.backgroundMesh = obj.mesh;
            sUm.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(sUm);
            uMesh.compute(levelSet);
        end

        function uMesh = computeHinge2Domain(obj,r)
            s.type        = 'Circle';
            s.radius      = r*obj.mesh.computeMeanCellSize();
            s.xCoorCenter = 0.75;
            s.yCoorCenter = 0.5;
            g             = GeometricalFunction(s);
            lsFun         = g.computeLevelSetFunction(obj.mesh);
            levelSet      = lsFun.fValues;
            sUm.backgroundMesh = obj.mesh;
            sUm.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(sUm);
            uMesh.compute(levelSet);
        end

        function uMesh = computeHinge3Domain(obj,r)
            s.type        = 'Circle';
            s.radius      = r*obj.mesh.computeMeanCellSize();
            s.xCoorCenter = 0.5;
            s.yCoorCenter = 0.25;
            g             = GeometricalFunction(s);
            lsFun         = g.computeLevelSetFunction(obj.mesh);
            levelSet      = lsFun.fValues;
            sUm.backgroundMesh = obj.mesh;
            sUm.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(sUm);
            uMesh.compute(levelSet);
        end

        function uMesh = computeHinge4Domain(obj,r)
            s.type        = 'Circle';
            s.radius      = r*obj.mesh.computeMeanCellSize();
            s.xCoorCenter = 0.25;
            s.yCoorCenter = 0.5;
            g             = GeometricalFunction(s);
            lsFun         = g.computeLevelSetFunction(obj.mesh);
            levelSet      = lsFun.fValues;
            sUm.backgroundMesh = obj.mesh;
            sUm.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(sUm);
            uMesh.compute(levelSet);
        end

        function filterPerimeter = createFilterPerimeter(obj)
            s.filterType    = 'PDE';
            s.boundaryType  = 'Neumann';
            s.mesh          = obj.mesh;
            s.trial         = LagrangianFunction.create(obj.mesh,1,'P1');
            f               = Filter.create(s);
            filterPerimeter = f;
        end

        function createGlobalPerimeterConstraint(obj)
            s.mesh              = obj.mesh;
            s.filter            = createFilterPerimeter(obj);
            s.epsilon           = 6*obj.mesh.computeMeanCellSize();
            s.value0            = 6;
            obj.globalPerimeter = PerimeterFunctional(s);
        end

        function createCost(obj,w)
            s.shapeFunctions{1} = obj.compliance;
            s.shapeFunctions{2} = obj.globalPerimeter;
            s.weights           = [1,w];
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            s.test  = LagrangianFunction.create(obj.mesh,1,'P1');
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh  = obj.mesh;
            s.type  = 'MassMatrix';
            LHS = LHSIntegrator.create(s);
            M = LHS.compute;

            h = obj.mesh.computeMinCellSize();
            M = h^2*eye(size(M));
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.shapeFunctions{2} = obj.hinge1Perimeter;
            s.shapeFunctions{3} = obj.hinge2Perimeter;
            s.shapeFunctions{4} = obj.hinge3Perimeter;
            s.shapeFunctions{5} = obj.hinge4Perimeter;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createDualVariable(obj)
            s.nConstraints   = 5;
            l                = DualVariable(s);
            obj.dualVariable = l;
        end

        function createOptimizer(obj,gJ)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 2000;
            s.tolerance      = 1e-8;
            s.constraintCase = {'INEQUALITY','INEQUALITY','INEQUALITY','INEQUALITY','INEQUALITY'};
            s.primal         = 'SLERP';                  
            s.etaNorm        = 0.02;
            s.etaNormMin     = 0.002;
            s.gJFlowRatio    = gJ;
            s.etaMax         = 0.1;
            s.etaMaxMin      = 0.01;
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