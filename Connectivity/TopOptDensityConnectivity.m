classdef TopOptDensityConnectivity < handle

    properties (Access = private)
        mesh
        filter
        filterComp
        filterPerimeter
        filterAdjointComp
        filterConnect
        filterAdjointConnect
        designVariable
        materialInterpolator
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
        gJ
        lambda1min
        primalUpdater
        p
    end

    methods (Access = public)

        function obj = TopOptDensityConnectivity() 
            for p = [0.0,3.0,10.0] 
                for lambda1min = [0.6,0.8,1.0,2.0] 
                    obj.p = p;
                    obj.lambda1min = lambda1min;
                    obj.init()
                    obj.createMesh();
                    obj.createDesignVariable();
                    obj.createFilterPerimeter();
                    obj.createFilterCompliance();
                    obj.createFilterConnectivity();
                    obj.createMaterialInterpolator();
                    obj.createCompliance();
                    obj.createConductivityInterpolator();
                    obj.createMassInterpolator();
                    obj.createElasticProblem();
                    obj.createComplianceFromConstitutive();
                    obj.createCompliance();
                    obj.createPerimeter();
                    obj.createEigenValueConstraint();                             
                    obj.createVolumeConstraint();
                    obj.createCost();
                    obj.createConstraint();
                    obj.createPrimalUpdater();
                    obj.createOptimizer();
                end
            end
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            x1      = linspace(0,2.0,100);
            x2      = linspace(0,1.0,50);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh.create(s);
        end
    
        function createDesignVariable(obj)
            s.fHandle = @(x) 1.0*ones(size(x(1,:,:)));
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
            f.updateEpsilon(1*obj.mesh.computeMinCellSize())
            obj.filterPerimeter = f;           
        end

        function createFilterConnectivity(obj)
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filterConnect = f; 
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
            s.mesh                 = obj.mesh;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            m = Material.create(s);
        end

        function createConductivityInterpolator(obj) 
            s.interpolation  = 'SimpAllThermal';
            s.f0   = 1e-3;                                             
            s.f1   = 1;  
            s.dim  = '2D';
            a = MaterialInterpolator.create(s);
            obj.conductivityInterpolator = a;            
        end 

        function createMassInterpolator(obj)
            s.interpolation  = 'SIMPThermal';                              
            s.f0   = 1e-3;
            s.f1   = 1;
            s.pExp = 1;
            a = MaterialInterpolator.create(s);
            obj.massInterpolator = a;            
        end      

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.boundaryConditions = obj.createElasticBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function c = createComplianceFromConstitutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstitutiveTensor(s);
        end

        function createCompliance(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filterComp;
            if ~isempty(obj.filterAdjointComp)
                s.filterAdjoint               = obj.filterAdjointComp;
            end
            s.complainceFromConstitutive  = obj.createComplianceFromConstitutive();
            s.material                    = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function createPerimeter(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filterComp;
            p = SimplePerimeterFunctional(s);
            obj.perimeter = p;
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
            s.targetEigenValue  = obj.lambda1min;    
            s.isCompl           = true;
            obj.minimumEigenValue = StiffnessEigenModesConstraint(s);
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.shapeFunctions{2} = obj.minimumEigenValue; 
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.shapeFunctions{2} = obj.perimeter;
%             s.shapeFunctions{3} = obj.minimumEigenValue;
            s.weights           = [1.0,obj.p];
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            test  = LagrangianFunction.create(obj.mesh,1,'P1');
            trial = LagrangianFunction.create(obj.mesh,1,'P1');
            M = IntegrateLHS(@(u,v) DP(v,u), test, trial, obj.mesh, 'Domain', 2);
            h = obj.mesh.computeMinCellSize();
            M = h^2*eye(size(M));
        end

        function createPrimalUpdater(obj)
            s.ub     = 1;
            s.lb     = 0;
            s.tauMax = 1000;
            s.tau    = [];
            obj.primalUpdater = ProjectedGradient(s);
        end

        function createOptimizer(obj)
%             s.monitoring     = true;
%             s.cost           = obj.cost;
%             s.constraint     = obj.constraint;
%             s.designVariable = obj.designVariable;
%             s.dualVariable   = obj.dualVariable;
%             s.maxIter        = 1000;
%             s.tolerance      = 1e-8;
%             s.constraintCase{1} = 'EQUALITY';
%             s.constraintCase{2} = 'INEQUALITY';                             
%             s.ub             = 1;
%             s.lb             = 0;
%             opt              = OptimizerMMA(s);
%             opt.solveProblem();
%             obj.optimizer = opt;
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.GIFname        = 'name';
            s.designVariable = obj.designVariable;
            s.maxIter        = 2000;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY','INEQUALITY'};
            s.primal         = 'PROJECTED GRADIENT';
            s.etaNorm        = 0.01;
            s.gJFlowRatio    = 2;
            s.primalUpdater  = obj.primalUpdater;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;

%             s.primal         = 'SLERP';
%             s.ub             = inf;
%             s.lb             = -inf;
%             s.etaNorm        = 0.02;
%             s.etaNormMin     = 0.02;
%             s.gJFlowRatio    = obj.gJ;
%             s.etaMax         = 1;
%             s.etaMaxMin      = 0.01;
%             s.filter         = obj.filterComp;
%             opt = OptimizerNullSpace(s);
%             opt.solveProblem();
%             obj.optimizer = opt;
            saveas(figure(1),'_p_'+string(obj.p)+string(obj.lambda1min)+'design.png','png')
            saveas(figure(2),'_p_'+string(obj.p)+string(obj.lambda1min)+'graficos.png','png')
            writematrix(obj.designVariable.fun.fValues,'_p_'+string(obj.p)+string(obj.lambda1min)+'.txt')

        end

        function bc = createElasticBoundaryConditions(obj)
            type = 'cantilever';
            if isequal(type, 'cantilever')
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
            elseif isequal(type, 'bridge')
                xMax    = max(obj.mesh.coord(:,1));
                yMax    = max(obj.mesh.coord(:,2));
                isDir   = @(coor)  abs(coor(:,2))==0.0 & (abs(coor(:,1))>= 0.95*xMax | abs(coor(:,1))<= 0.05*xMax);
                isForce = @(coor)  (abs(coor(:,2))==yMax & abs(coor(:,1))>=0.45*xMax & abs(coor(:,1))<=0.55*xMax);
    
                sDir{1}.domain    = @(coor) isDir(coor);
                sDir{1}.direction = [1,2];
                sDir{1}.value     = 0;
    
                sPL{1}.domain    = @(coor) isForce(coor);
                sPL{1}.direction = 2;
                sPL{1}.value     = -1;
            elseif isequal(type, 'acantilever')
                xMax    = max(obj.mesh.coord(:,1));
                yMax    = max(obj.mesh.coord(:,2));
                isDir   = @(coor)  abs(coor(:,1))== 0.0;
                isForce = @(coor)  abs(coor(:,1)) == xMax & abs(coor(:,2))==0.0;
    
                sDir{1}.domain    = @(coor) isDir(coor);
                sDir{1}.direction = [1,2];
                sDir{1}.value     = 0;
    
                sPL{1}.domain    = @(coor) isForce(coor);
                sPL{1}.direction = 2;
                sPL{1}.value     = -1;
            end
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