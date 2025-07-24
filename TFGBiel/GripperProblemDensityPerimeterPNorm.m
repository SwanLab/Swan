classdef GripperProblemDensityPerimeterPNorm < handle

    properties (Access = private)
        filename
        mesh
        filter
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        volume
        perimeter
        globalPerimeter
        cost
        constraint
        dualVariable
        optimizer
    end

    methods (Access = public)

        function obj = GripperProblemDensityPerimeterPNorm(p,pT,gJ,w)
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createNonSelfAdjCompliance();
            obj.createVolumeConstraint();
            obj.createPerimeterConstraint(p,pT);
            obj.createGlobalPerimeterConstraint();
            obj.createCost(w);
            obj.createConstraint();
            obj.createDualVariable();
            obj.createOptimizer(gJ);

            fileLocation = 'C:\Users\Biel\Desktop\UNI\TFG\ResultatsNormP_Density\00. From Batch';
            
            vtuName = fullfile(fileLocation, sprintf('Topology_Gripper_perimeter_p%d_ptarget%.2f_gJ%.2f_eta0.02_w%.2f',p,pT,gJ,w));
            obj.designVariable.fun.print(vtuName);
            
            figure(2)
            set(gcf, 'Position', get(0, 'Screensize'));
            fileName1 = fullfile(fileLocation, sprintf('Monitoring_Gripper_perimeter_p%d_ptarget%.2f_gJ%.2f_eta0.02_w%.2f.fig',p,pT,gJ,w));
            fileName2 = fullfile(fileLocation, sprintf('Monitoring_Gripper_perimeter_p%d_ptarget%.2f_gJ%.2f_eta0.02_w%.2f.png',p,pT,gJ,w));
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

        function createPerimeterConstraint(obj,p,pTarget)
            s.mesh            = obj.mesh;
            s.perimeterTarget = pTarget;
            s.p               = p;
            s.gradientTest    = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.perimeter     = PerimeterNormPFunctional(s);
        end

        function createGlobalPerimeterConstraint(obj)
            s.mesh              = obj.mesh;
            s.filter            = createFilterPerimeter(obj);
            s.epsilon           = 6*obj.mesh.computeMeanCellSize();
            s.value0            = 6;
            obj.globalPerimeter = PerimeterFunctional(s);
        end

        function filterPerimeter = createFilterPerimeter(obj)
            s.filterType    = 'PDE';
            s.boundaryType  = 'Robin';
            s.mesh          = obj.mesh;
            s.trial         = LagrangianFunction.create(obj.mesh,1,'P1');
            f               = Filter.create(s);
            filterPerimeter = f;
        end

        function createCost(obj,w)
            s.shapeFunctions{1} = obj.compliance;
            s.shapeFunctions{2} = obj.globalPerimeter;
            s.weights           = [1,w];
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
            s.shapeFunctions{2} = obj.perimeter;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createDualVariable(obj)
            s.nConstraints   = 2;
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
            s.constraintCase = {'INEQUALITY','INEQUALITY'};
            s.primal         = 'PROJECTED GRADIENT';
            s.tauMax         = Inf;
            s.ub             = 1;
            s.lb             = 0;
            s.etaNorm        = 0.02;
            s.gJFlowRatio    = gJ;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function m = createMaterial(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = obj.filter.compute(f{1},1);            
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            s.mesh = obj.mesh;
            m = Material.create(s);
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