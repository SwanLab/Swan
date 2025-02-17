classdef TopOptDensityConnectivity < handle

    properties (Access = private)
        mesh
        filter
        filterComp
        filterAdjointComp
        filterConnect
        filterAdjointConnect
        designVariable
        materialInterpolator
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
    end

    methods (Access = public)

        function obj = TopOptDensityConnectivity()
%             for gJ = [0.2, 1.0, 2.0]
                for lambda1min = [0.05]
%                     obj.gJ = gJ;
                    obj.lambda1min = lambda1min;        
                    obj.init()
                    obj.createMesh();
                    obj.createDesignVariable();
                    obj.createFilterCompliance();
                    obj.createFilterConnectivity();
                    obj.createMaterialInterpolator();
                    obj.createElasticProblem();
                    obj.createComplianceFromConstitutive();
                    obj.createCompliance();
                    obj.createPerimeter();
                    obj.createEigenValueConstraint();                             
                    obj.createVolumeConstraint();
                    obj.createCost();
                    obj.createConstraint();
                    obj.createDualVariable();
                    obj.createOptimizer();
                end
%             end
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
%             s.filterType = 'FilterAndProject';
%             s.mesh       = obj.mesh;
%             s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.filterStep = 'LUMP';
%             s.beta       = 2.0;
% %             s.eta        = 0.5;
%             obj.filterComp = Filter.create(s);
%             
%             s.filterType = 'FilterAdjointAndProject';   
%             s.mesh       = obj.mesh;
%             s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.filterStep = 'LUMP';
%             s.beta       = 2.0;
% %             s.eta        = 0.5;
%             obj.filterAdjointComp = Filter.create(s);

            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filterComp = f;
        end
        

        function createFilterConnectivity(obj)
%             s.filterType = 'FilterAndProject';
%             s.mesh       = obj.mesh;
%             s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.filterStep = 'LUMP';
%             s.beta       = 6.0;
%             s.eta        = 0.5;
%             obj.filterConnect = Filter.create(s);
% 
%             s.filterType = 'FilterAdjointAndProject';   
%             s.mesh       = obj.mesh;
%             s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.filterStep = 'LUMP';
%             s.beta       = 6.0;
%             s.eta        = 0.5;
%             obj.filterAdjointConnect = Filter.create(s);

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

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filterComp;
            if ~isempty(obj.filterAdjointComp)
                s.filterAdjoint               = obj.filterAdjointComp;
            end
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createEigenValueConstraint(obj)                           
            s.mesh              = obj.mesh;
            s.designVariable    = obj.designVariable;
            s.filter            = obj.filterConnect;
            s.filterAdjoint     = obj.filterAdjointConnect;   
            s.targetEigenValue  = obj.lambda1min;      
            s.shift             = 0.0;
            obj.minimumEigenValue = StiffnesEigenModesConstraint(s);
%             eigen = StiffnessEigenModesComputer(s);
%             s.eigenModes = eigen;
%             obj.minimumEigenValue = MinimumEigenValueFunctional(s);
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
            s.weights           = [1.0; 0.0]; %; -1.0];
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
        end

        function createDualVariable(obj)
            s.nConstraints   = 2;                                        
            l                = DualVariable(s);
            obj.dualVariable = l;
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 1000;
            s.tolerance      = 1e-8;
            s.constraintCase{1} = 'EQUALITY';
            s.constraintCase{2} = 'INEQUALITY';                             
            s.ub             = 1;
            s.lb             = 0;
            opt              = OptimizerMMA(s);
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
            saveas(figure(1),'DEN1e-3lambda1min'+string(obj.lambda1min)+'design.png','png')
            saveas(figure(2),'DEN1e-3lambda1min'+string(obj.lambda1min)+'graficos.png','png')
            writematrix(obj.designVariable.fun.fValues,'DEN1e-3lambda1min'+string(obj.lambda1min)+'.txt')
        end

        function bc = createBoundaryConditions(obj)
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
    end
end
