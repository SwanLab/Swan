classdef TopOpt_Density_Chair< handle

    properties (Access = private)
        mesh
        rbe3Mesh
        rbe3Nodes
        filename
        filter
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        volume
        cost
        constraint
        primalUpdater
        optimizer
    end

    methods (Access = public)

        function obj = TopOpt_Density_Chair(vTar)
            filename = 'Chair';
            obj.init()
            obj.createMesh(filename);
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createVolumeConstraint(vTar);
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();

            d = obj.designVariable.fun;

            saveas(gcf,['ChairSupport/Results/Vol',num2str(vTar),'/MonitoringDensity_',filename,'.fig']);
            d.print(['ChairSupport/Results/Vol',num2str(vTar),'/Density',filename,'_fValues']);

            fem = obj.physicalProblem;
            sigma = fem.stressFun;
            devSig = Deviatoric(sigma);
            vonMises = sqrt(1.5.*DDP(devSig,devSig));

            sF.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            sF.mesh = obj.mesh;
            sF.filterType = 'PDE';
            filter = Filter.create(sF);
            filter.updateEpsilon(2*obj.mesh.computeMeanCellSize());
            vmSig = filter.compute(vonMises,3);
            vmSig.print(['ChairSupport/Results/Vol',num2str(vTar),'/Density_',filename,'_VMFinal'])
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj,fName)
            file = fName;
            obj.filename = file;
            a.fileName = file;
            s = FemDataContainer(a);
            s.mesh.readBoundaryMeshFromGiD([fName,'_Boundary'],0);
            obj.mesh = s.mesh;
        end

        function createDesignVariable(obj)
            s.fHandle  = @(x) ones(size(x(1,:,:)));
            s.ndimf    = 1;
            s.mesh     = obj.mesh;
            aFun       = AnalyticalFunction(s);
            s.fun      = aFun.project('P1');
            s.type     = 'Density';
            s.plotting = true;

            % FIXED REGIONS
            region1 = @(x) (x(:,2) >= ...
                107.61 + ((58.471-107.61)/(-216.13+196.28))*(x(:,1)+196.28)) ...
                & (x(:,2) < 106.79) ...
                & (x(:,2) > 60.5328);

            region2 = @(x) x(:,2) <= ...
                72.7349 +  (40.492 - 72.7349)/(212.701 - 220.016)*(x(:,1) - 220.016)...
                & (x(:,2) < 72.4) ...
                & (x(:,2) > 38);

            region3 = @(x) x(:,2) <= ...
                140.94 + (110.23 - 140.94)/(238.43 - 251.88)*(x(:,1) - 251.88)...
                & (x(:,2) < 140.94) ...
                & (x(:,2) > 108.76);

            s.isFixed  = obj.computeFixedVolumeDomain(...
                         @(x) region1(x) | region2(x)| region3(x));
            
            ls     = DesignVariable.create(s);
            obj.designVariable = ls;
        end

        function isFixed = computeFixedVolumeDomain(obj,cond)
            coor    = obj.mesh.coord;
            isFixed = find(cond(coor));
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

        function createMaterialInterpolator(obj)
            E0   = 7;
            nu0  = 1/3;
            E1   = 70e3;
            nu1  = 1/3;
            ndim =  obj.mesh.ndim;

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
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filter;
            s.complainceFromConstitutive = obj.createComplianceFromConstiutive();
            s.material                   = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function uMesh = createBaseDomain(obj)
            sG.type          = 'Full';
            g                = GeometricalFunction(sG);
            lsFun            = g.computeLevelSetFunction(obj.mesh);
            levelSet         = lsFun.fValues;
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh            = UnfittedMesh(s);
            uMesh.compute(levelSet);
        end

        function createVolumeConstraint(obj,vTar)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.test = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = vTar;
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
            M = IntegrateLHS(@(u,v) u.*v,obj.designVariable.fun,obj.designVariable.fun,obj.mesh,'Domain',3);
            M = diag(sum(M,1));
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createPrimalUpdater(obj)
            rho      = obj.designVariable;
            fixedDof = rho.getFixedDofs();
            s.ub     = ones(size(rho.fun.fValues));
            s.lb     = zeros(size(rho.fun.fValues));
            s.lb(fixedDof) = 1;
            s.tauMax = 500;
            s.tau    = [];
            obj.primalUpdater = ProjectedGradient(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 500;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.primal         = 'PROJECTED GRADIENT';
            s.primalUpdater  = obj.primalUpdater;
            s.delta          = 0.005;
            s.deltaMin       = 0.001;
            s.etaStar        = 0.6;
            s.etaMax0        = 0.1;
            s.etaMaxMin      = 0.0005;
            s.gif            = false;
            s.gifName        = 'CadiraGif';
            s.printing       = false;
            s.printName      = [];
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function m = createMaterial(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = f{1}.project('P1');
            % f = obj.filter.compute(f{1},1);            
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            s.mesh                 = obj.mesh;
            m = Material.create(s);
        end

        function bc = createBoundaryConditions(obj)
            femReader = FemInputReaderGiD();
            s         = femReader.read(obj.filename);
            sPL       = obj.computeCondition(s.pointload);
            sDir     = obj.computeCondition(s.dirichlet);

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