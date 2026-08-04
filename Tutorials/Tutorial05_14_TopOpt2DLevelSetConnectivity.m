classdef Tutorial05_14_TopOpt2DLevelSetConnectivity< handle

    properties (Access = private)
        mesh
        filterPerimeter
        filterComp
        filterConnect
        filterAdjointConnect
        designVariable
        C
        dC
        conductivityInterpolator
        massInterpolator
        physicalProblem
        compliance
        volume
        cost
        constraint
        dualVariable
        optimizer
        minimumEigenValue
        perimeter
        primalUpdater
        E0; nu0; E1; nu1
    end 

    methods (Access = public)
        function obj = Tutorial05_14_TopOpt2DLevelSetConnectivity()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilterPerimeter();
            obj.createFilterCompliance();
            obj.createFilterConnectivity();
            obj.createMaterial();
            obj.createElasticProblem();
            obj.createComplianceFromConstitutive();
            obj.createCompliance();
            obj.createConductivityInterpolator();
            obj.createMassInterpolator();
            obj.createEigenValueConstraint();   
            obj.createPerimeter();                  
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
            obj.E0  = 1e-5;
            obj.nu0 = 1/3;
            obj.E1  = 1;
            obj.nu1 = 1/3;
        end

        function createMesh(obj)
            obj.mesh = TriangleMesh(2,1,100,50);
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

        function createFilterCompliance(obj)
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filterComp = f;           
        end

       function createFilterPerimeter(obj)
            s.filterType = 'PDE';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filterPerimeter = f;           
        end
 
        function createFilterConnectivity(obj)
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filterConnect = f;
            obj.filterAdjointConnect =[];
        end

        function createConductivityInterpolator(obj) 
            s.f0   = 1e-3;                                             
            s.f1   = 1;  
            s.dim  = '2D';
            a = SimpAllThermalInterpolation(s);
            obj.conductivityInterpolator = a;            
        end 

        function createMassInterpolator(obj)
            s.f0   = 1e-3;
            s.f1   = 1;
            s.pExp = 1;
            a = SIMPThermalInterpolation(s);
            obj.massInterpolator = a;            
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

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = [];
            s.dim = '2D';
            s.boundaryConditions = obj.createElasticBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = DirectSolver();
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function c = createComplianceFromConstitutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstitutiveTensor(s);
        end

        function createCompliance(obj)
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filterComp;
            s.complainceFromConstitutive = obj.createComplianceFromConstitutive();
            s.C                          = obj.C;
            s.dC                         = obj.dC;
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


        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.test = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;
            s.uMesh = obj.createBaseDomain();
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createEigenValueConstraint(obj)                           
            s.mesh              = obj.mesh;
            s.designVariable    = obj.designVariable;
            s.filter            = obj.filterConnect;
            s.filterAdjoint     = obj.filterAdjointConnect;
            s.boundaryConditions = obj.createEigenvalueBoundaryConditions();
            s.conductivityInterpolator = obj.conductivityInterpolator; 
            s.massInterpolator         = obj.massInterpolator; 
            s.targetEigenValue  = 1.0;    
            s.isComplementary   = true; % 1 - chi
            obj.minimumEigenValue = StiffnessEigenModesConstraint(s);
        end

        function createPerimeter(obj)
            s.mesh        = obj.mesh;
            s.filter      = obj.filterPerimeter;
            s.epsilon     = obj.mesh.computeMeanCellSize();
            s.value0      = 6;
            s.uMesh       = obj.createBaseDomain();
            P             = PerimeterFunctional(s);
            obj.perimeter = P;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.shapeFunctions{2} = obj.perimeter;
            s.weights           = [1.0,0.5]; 
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
            s.shapeFunctions{2} = obj.minimumEigenValue; 
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createPrimalUpdater(obj)
            s.mesh = obj.mesh;
            obj.primalUpdater = SLERP(s);
        end

        function createOptimizer(obj,max)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 3;
            s.tolerance      = 1e-8;
            s.constraintCase{1} = 'EQUALITY';
            s.constraintCase{2} = 'INEQUALITY';
            s.primalUpdater  = obj.primalUpdater;
            s.delta          = 0.02;
            s.deltaMin       = 0.02;
            s.etaStar        = 0.2;
            s.etaMax0        = 1;
            s.etaMaxMin      = 0.01;
            s.gif            = false;
            s.gifName        = 'Tutorial05_3';
            s.printing       = false;
            s.printName      = 'Tutorial05_3';
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createElasticBoundaryConditions(obj)
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

        function  bc = createEigenvalueBoundaryConditions(obj)
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isLeft  = @(coor) abs(coor(:,1))==xMin;
            isRight = @(coor) abs(coor(:,1))==xMax;
            isFront = @(coor) abs(coor(:,2))==yMin;
            isBack = @(coor) abs(coor(:,2))== yMax;
            isDir   = @(coor) isLeft(coor) | isRight(coor) | isFront(coor) | isBack(coor);  
            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = 1;
            sDir{1}.value     = 0;
            sDir{1}.ndim = 1;
            
            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;
            s.pointloadFun = [];

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);  
        end

    end
end