classdef TopOpt_LevelSet_Puente < handle

    properties (Access = private)
        mesh
        filename
        filter
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
        E0; nu0; E1; nu1
    end

    methods (Access = public)

        function obj = TopOpt_LevelSet_Puente(vTar)
            filename = 'Puente';
            obj.init()
            obj.createMesh(filename);
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterial();
            obj.createElasticProblem(filename,vTar);
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createVolumeConstraint(vTar);
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();

            d = obj.designVariable.fun;

            saveas(gcf,['CIM_UPC/Results/Vol',num2str(vTar),'/MonitoringLevelSet_',filename,'.fig']);
            d.print(['CIM_UPC/Results/Vol',num2str(vTar),'/LevelSet',filename,'_fValues']);

            fem = obj.physicalProblem;
            uF = fem.uFun;
            uF.print(['CIM_UPC/Results/Vol',num2str(vTar),'/LevelSet_',fName,'_DisplFinal']);
            sigma = fem.stressFun;
            devSig = Deviatoric(sigma);
            vonMises = sqrt(1.5.*DDP(devSig,devSig));

            sF.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            sF.mesh = obj.mesh;
            sF.filterType = 'PDE';
            filter = Filter.create(sF);
            filter.updateEpsilon(3*obj.mesh.computeMeanCellSize());
            vmSig = filter.compute(vonMises,3);
            vmSig.print(['CIM_UPC/Results/Vol',num2str(vTar),'/LevelSet_',filename,'_VMFinal'])
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
            obj.E0  = 70;
            obj.nu0 = 1/3;
            obj.E1  = 70e3;
            obj.nu1 = 1/3;
        end

        function createMesh(obj,fName)
            file = fName;
            obj.filename = file;
            a.fileName = file;
            s = FemDataContainer(a);
            s.mesh.readBoundaryMeshFromGiD([fName,'_Boundary'],-1);
            obj.mesh = s.mesh;
        end

        function createDesignVariable(obj)
            s.type = 'Full';
            g      = GeometricalFunction(s);
            lsFun  = g.computeLevelSetFunction(obj.mesh);
            s.fun  = lsFun;
            s.mesh = obj.mesh;
            s.type = 'LevelSet';
            s.plotting = false;
            ls     = DesignVariable.create(s);
            obj.designVariable = ls;
        end

        function createFilter(obj)
            s.filterType = 'PDE';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            f.updateEpsilon(3*obj.mesh.computeMeanCellSize());
            obj.filter = f;
        end

        function createMaterial(obj)
            N = obj.mesh.ndim;
            muA    = obj.computeMu(obj.E0,obj.nu0);
            kappaA = obj.computeKappa(obj.E0,obj.nu0,N);
            muB    = obj.computeMu(obj.E1,obj.nu1);
            kappaB = obj.computeKappa(obj.E1,obj.nu1,N);

            mu    = @(rho) SimpAllInterpolator.computeMu(muA,muB,kappaA,kappaB,rho,N); mu = @(rho) Expand(mu(rho),4);
            kappa  = @(rho) SimpAllInterpolator.computeKappa(muA,muB,kappaA,kappaB,rho,N); kappa = @(rho) Expand(kappa(rho),4);
            lambda = @(rho) kappa(rho) - (2/N)*mu(rho);
            I      = ConstantFunction.create(eye4D(N),obj.mesh);
            IxI    = ConstantFunction.create(kronEye(N),obj.mesh);
            obj.C  = @(rho) 2*mu(rho{1}).*I + lambda(rho{1}).*IxI;

            dmu     = @(rho) SimpAllInterpolator.computeMuDerivative(muA,muB,kappaA,kappaB,rho,N); dmu = @(rho) Expand(dmu(rho),4);
            dkappa  = @(rho) SimpAllInterpolator.computeKappaDerivative(muA,muB,kappaA,kappaB,rho,N); dkappa = @(rho) Expand(dkappa(rho),4);
            dlambda = @(rho) dkappa(rho) - (2/N)*dmu(rho);
            obj.dC  = @(rho) {2*dmu(rho{1}).*I + dlambda(rho{1}).*IxI};
         end

        function mu = computeMu(obj,E,nu)
            mu = E./(2*(1+nu));
        end

        function kappa = computeKappa(obj,E,nu,N)
            kappa = E./(N*(1-(N-1)*nu));
        end

        function createElasticProblem(obj,fName,vTar)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = [];
            s.dim = '3D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = DirectSolver();
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;

            N = obj.mesh.ndim;
            E = obj.E1; nu = obj.nu1; lam = LameParametersConverter;
            mu     = lam.computeShearFromYoungAndPoisson(E,nu);  mu = Expand(mu,4);
            lambda = lam.computeLambdaFromYoungAndPoisson(E,nu,N);  lambda = Expand(lambda,4);
            I      = ConstantFunction.create(eye4D(N),obj.mesh);
            IxI    = ConstantFunction.create(kronEye(N),obj.mesh);
            m = 2*mu.*I + lambda.*IxI;
            fem.updateMaterial(m);
            fem.solve();

            uF = fem.uFun;
            uF.print(['CIM_UPC/Results/Vol',num2str(vTar),'/LevelSet_',fName,'_DisplInitial']);

            sigma = fem.stressFun;
            devSig = Deviatoric(sigma);
            vonMises = sqrt(1.5.*DDP(devSig,devSig));

            sF.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            sF.mesh = obj.mesh;
            sF.filterType = 'PDE';
            filt = Filter.create(sF);
            filt.updateEpsilon(3*obj.mesh.computeMeanCellSize());
            vmSig = filt.compute(vonMises,3);
            vmSig.print(['CIM_UPC/Results/Vol',num2str(vTar),'/LevelSet_',fName,'_VMInitial']);
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
            s.mesh = obj.mesh;
            obj.primalUpdater = SLERP(s);
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
            s.delta          = 0.001;
            s.deltaMin       = 0.001;
            s.etaStar        = 1.0;
            s.etaMax0        = 0.5;
            s.etaMaxMin      = 0.01;
            s.gif            = false;
            s.gifName        = [];
            s.printing       = false;
            s.printName      = [];
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj)
            femReader = FemInputReaderGiD();
            s         = femReader.read(obj.filename);
            s.pointload(:,3) = 300*s.pointload(:,3); % To capture 6.4mm "flecha"
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