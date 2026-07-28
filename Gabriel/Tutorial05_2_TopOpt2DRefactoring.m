classdef Tutorial05_2_TopOpt2DRefactoring < handle

    properties (Access = private)
        mesh
        filter
        designVariable
        material
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

        function obj = Tutorial05_2_TopOpt2DRefactoring()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterial();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
            obj.E0  = 1e-3;
            obj.nu0 = 1/3;
            obj.E1  = 1;
            obj.nu1 = 1/3;
        end

        function createMesh(obj)
            obj.mesh = TriangleMesh(2,1,100,50);
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);

            sD.fun      = aFun.project('P1');
            sD.mesh     = obj.mesh;
            sD.type     = 'Density';
            sD.plotting = true;
            dens        = DesignVariable.create(sD);
            obj.designVariable = dens;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

        function createMaterial(obj)
            s.tensorFcn           = @(rho) obj.computeTensor(rho);
            s.tensorDerivativeFcn = @(rho) obj.computeTensorDerivative(rho);
            obj.material = MaterialFromFunctions(s);
        end

       

        function mu = computeMu(obj,E,nu)
            mu = E./(2*(1+nu));
        end

        function kappa = computeKappa(obj,E,nu,N)
            kappa = E./(N*(1-(N-1)*nu));
        end

      

        function etaMu = computeEtaMu(obj,mu,kappa,N)
            num = -mu.*(4.*mu - kappa.*N^2 - 2.*mu.*N^2 + 2.*mu.*N);
            den = 2.*N.*(kappa + 2*mu);
            etaMu = num./den;
        end

        function etaKappa = computeEtaKappa(obj,mu,N)
            etaKappa = 2.*mu.*(N-1)./N;
        end

        function c = computeCoeff(obj,f0,f1,eta0,eta1)
            c.n01 = -(f1 - f0).*(eta1 - eta0);
            c.n0  = f0.*(f1 + eta0);
            c.n1  = f1.*(f0 + eta1);
            c.d0  = (f1 + eta0);
            c.d1  = (f0 + eta1);
        end

        function f = computeRationalFunction(obj,rho,c)
            num = c.n01.*(1-rho).*rho + c.n0.*(1-rho) + c.n1.*rho;
            den = c.d0.*(1-rho) + c.d1.*rho;
            f   = num./den;
        end

        function df = computeRationalDerivative(obj,rho,c)
            num    = c.n01.*(1-rho).*rho + c.n0.*(1-rho) + c.n1.*rho;
            den    = c.d0.*(1-rho) + c.d1.*rho;
            derNum = (c.n01.*(1-2.*rho)-c.n0+c.n1).*den - num.*(-c.d0+c.d1);
            df     = derNum./(den.^2);
        end

        function [mu,kappa] = computeMuKappa(obj,rho)
            N = obj.mesh.ndim;
            m0 = obj.computeMu(obj.E0,obj.nu0);   k0 = obj.computeKappa(obj.E0,obj.nu0,N);
            m1 = obj.computeMu(obj.E1,obj.nu1);   k1 = obj.computeKappa(obj.E1,obj.nu1,N);

            eta0mu = obj.computeEtaMu(m0,k0,N);
            eta1mu = obj.computeEtaMu(m1,k1,N);
            cMu    = obj.computeCoeff(m0,m1,eta0mu,eta1mu);
            mu     = obj.computeRationalFunction(rho,cMu);

            eta0k  = obj.computeEtaKappa(m0,N);
            eta1k  = obj.computeEtaKappa(m1,N);
            cKappa = obj.computeCoeff(k0,k1,eta0k,eta1k);
            kappa  = obj.computeRationalFunction(rho,cKappa);
        end

        function [dmu,dkappa] = computeMuKappaDerivative(obj,rho)
            N = obj.mesh.ndim;
            m0 = obj.computeMu(obj.E0,obj.nu0);   k0 = obj.computeKappa(obj.E0,obj.nu0,N);
            m1 = obj.computeMu(obj.E1,obj.nu1);   k1 = obj.computeKappa(obj.E1,obj.nu1,N);

            eta0mu = obj.computeEtaMu(m0,k0,N);
            eta1mu = obj.computeEtaMu(m1,k1,N);
            cMu    = obj.computeCoeff(m0,m1,eta0mu,eta1mu);
            dmu    = obj.computeRationalDerivative(rho,cMu);

            eta0k  = obj.computeEtaKappa(m0,N);
            eta1k  = obj.computeEtaKappa(m1,N);
            cKappa = obj.computeCoeff(k0,k1,eta0k,eta1k);
            dkappa = obj.computeRationalDerivative(rho,cKappa);
        end

        

        function C = computeTensor(obj,rho)
            N = obj.mesh.ndim;
            [mu,kappa] = obj.computeMuKappa(rho);
            lambda = kappa - (2/N).*mu;
            mu     = Expand(mu,4);
            lambda = Expand(lambda,4);
            I   = ConstantFunction.create(eye4D(N),obj.mesh);
            IxI = ConstantFunction.create(kronEye(N),obj.mesh);
            C = (mu.*I).*2 + lambda.*IxI;  
        end

        function dC = computeTensorDerivative(obj,rho)
            N = obj.mesh.ndim;
            [dmu,dkappa] = obj.computeMuKappaDerivative(rho);
            dlambda = dkappa - (2/N).*dmu;
            dmu     = Expand(dmu,4);
            dlambda = Expand(dlambda,4);
            I   = ConstantFunction.create(eye4D(N),obj.mesh);
            IxI = ConstantFunction.create(kronEye(N),obj.mesh);
            dC = (dmu.*I).*2 + dlambda.*IxI;   
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.material;
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
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filter;
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.material;
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function uMesh = createBaseDomain(obj)
            levelSet         = -ones(obj.mesh.nnodes,1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
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

        function createPrimalUpdater(obj)
            s.ub     = 1;
            s.lb     = 0;
            s.tauMax = 1000;
            s.tau    = [];
            obj.primalUpdater = ProjectedGradient(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 3;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.primal         = 'PROJECTED GRADIENT';
            s.etaNorm        = 0.01;
            s.gJFlowRatio    = 2;
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
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.4*yMax & abs(coor(:,2))<=0.6*yMax);

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = -1;

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
end