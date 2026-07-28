classdef TopOpt_Density_C5Sup < handle

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

        function obj = TopOpt_Density_C5Sup(vTar)
            filename = 'C5_Sup';
            obj.init()
            obj.createMesh(filename);
            obj.createRBE3Meshes(filename);
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem(filename,vTar);
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createVolumeConstraint(vTar);
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();

            d = obj.project();

            saveas(gcf,['Addicraft/Results/Vol',num2str(vTar),'/MonitoringDensity_',filename,'.fig']);
            d.print(['Addicraft/Results/Vol',num2str(vTar),'/Density_',filename,'_fValues']);

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
            vmSig.print(['Addicraft/Results/Vol',num2str(vTar),'/Density_',filename,'_VMFinal'])
        end

        function d = project(obj)
            s.mesh = obj.mesh;
            s.trial = obj.designVariable.fun;
            s.filterStep = 'PDE';
            s.eta = 0.3;
            s.beta = 10;
            filt = FilterAndProject(s);
            d = filt.compute(obj.designVariable.fun,3);
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
            s.mesh.readBoundaryMeshFromGiD([fName,'_Boundary'],-1);
            obj.mesh = s.mesh;
        end

        function createRBE3Meshes(obj,fName)
            for i=1:3
                run(fName);
                rollerRowsEl = External_border_elements(:,1)==i-1;
                rollerRowsN  = External_border_nodes(:,1)==i-1;
                s.coord = obj.mesh.coord;
                s.connec = External_border_elements(rollerRowsEl,3:5);
                s.kFace = -1;
                s.type = 'TRIANGLE';
                bDirMesh = SurfaceMesh(s);
                obj.rbe3Mesh{i} = bDirMesh.computeCanonicalMesh();
                obj.rbe3Nodes{i} = External_border_nodes(rollerRowsN,2);
            end
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);

            fixN = [];
            for i=1:3
                fixN = [fixN;obj.rbe3Nodes{i}];
            end
            
            sD.fun      = aFun.project('P1');
            sD.mesh     = obj.mesh;
            sD.type     = 'Density';
            sD.isFixed  = unique(fixN);
            sD.plotting = false;
            dens        = DesignVariable.create(sD);
            obj.designVariable = dens;
        end

        function createFilter(obj)
            s.filterType = 'PDE';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            f.updateEpsilon(2*obj.mesh.computeMeanCellSize());
            obj.filter = f;
        end

        function createMaterialInterpolator(obj)
            E0   = 70;
            nu0  = 1/3;
            E1   = 70e3;
            nu1  = 1/3;
            ndim = 3;

            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPALL';
            s.dim            = '3D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;
        end

        function createElasticProblem(obj,fName,vTar)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '3D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.rbe3Mesh = obj.rbe3Mesh;
            s.rbe3Nodes = obj.rbe3Nodes;
            s.rbe3Value(1,1) = 0;
            s.rbe3Value(1,2) = 0;
            s.rbe3Value(1,3) = 0;
            s.rbe3Value(2,1) = 0;
            s.rbe3Value(2,2) = 0;
            s.rbe3Value(2,3) = 0;
            s.rbe3Value(3,1) = NaN;
            s.rbe3Value(3,2) = NaN;
            s.rbe3Value(3,3) = 50;
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'ROLLER';
            s.solverCase = DirectSolver();
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;

            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = 3;
            s.young   = ConstantFunction.create(70e3,obj.mesh);
            s.poisson    = ConstantFunction.create(1/3,obj.mesh);
            m = Material.create(s);
            fem.updateMaterial(m);
            fem.solve();

            sigma = fem.stressFun;
            devSig = Deviatoric(sigma);
            vonMises = sqrt(1.5.*DDP(devSig,devSig));

            sF.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            sF.mesh = obj.mesh;
            sF.filterType = 'PDE';
            filt = Filter.create(sF);
            filt.updateEpsilon(2*obj.mesh.computeMeanCellSize());
            vmSig = filt.compute(vonMises,3);
            vmSig.print(['Addicraft/Results/Vol',num2str(vTar),'/Density_',fName,'_VMInitial']);
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
            s.maxIter        = 1000;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.primalUpdater  = obj.primalUpdater;
            s.etaNorm        = 0.01;
            s.etaNormMin     = 0.02;
            s.gJFlowRatio    = 0.5;
            s.etaMax         = 1;
            s.etaMaxMin      = 0.01;
            s.gif            = false;
            s.gifName        = [];
            s.printing       = false;
            s.printName      = [];
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
            s.dim                  = '3D';
            s.mesh                 = obj.mesh;
            m = Material.create(s);
        end

        function bc = createBoundaryConditions(obj)
            femReader = FemInputReaderGiD();
            s         = femReader.read(obj.filename);
            sPL       = obj.computeCondition(s.pointload);

            dirichletFun = [];

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