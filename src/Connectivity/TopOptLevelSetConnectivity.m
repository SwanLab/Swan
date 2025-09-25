classdef TopOptLevelSetConnectivity< handle

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
        filterAdjointProj 
        filterProj
        complianceProj
        volumeProj
        eta
        beta
        primalUpdater
        eigenvalue
    end 

    methods (Access = public)
        function obj = TopOptLevelSetConnectivity()
                for lambda1min = [0.4] %0.15]0.05,0.15,0.5,
                    obj.gJ = 1.0;
                    obj.lambda1min = lambda1min;
                    obj.init()
                    obj.createMesh();
                    obj.createDesignVariable();
                    obj.createFilter();
                    obj.createFilterConnectivity();
                    obj.createMaterialInterpolator();
                    obj.createElasticProblem();
                    obj.createComplianceFromConstitutive();
                    obj.createCompliance();
                    obj.createEigenValueConstraint();   
                    obj.createEigenValue()          
                    obj.createPerimeter();                  
                    obj.createVolumeConstraint();
                    obj.createCost();
                    obj.createConstraint();
                    obj.createPrimalUpdater();
                    obj.createOptimizer();
                end
%                 end
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
%             UnitMesh better
%             x1      = linspace(0,6,180);
%             x2      = linspace(0,1,30);
            x1      = linspace(0,2,120);
            x2      = linspace(0,1,60);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh.create(s);
        end

        function createDesignVariable(obj)
            s.type = 'Full';
%             s.type = 'CircleInclusion';
%             s.radius = 0.1;
%             s.xCoorCenter = 3.0;
%             s.yCoorCenter = 0.5;   
            g      = GeometricalFunction(s);
            lsFun  = g.computeLevelSetFunction(obj.mesh);
%             lsFun = LagrangianFunction.create(obj.mesh,1,'P1');
%             lsFun.setFValues(importdata('teste.txt'))
            s.fun  = lsFun;
            s.mesh = obj.mesh;
            s.type = 'LevelSet';
            s.plotting = true;
            ls     = DesignVariable.create(s);
%             ls.fun.setFValues(importdata('teste.txt'))
            obj.designVariable = ls;
            ls.fun.print('teste2D','Paraview')
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
% % 
%             s.filterType = 'FilterAndProject';
%             s.mesh       = obj.mesh;
%             s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.filterStep = 'PDE';
%             s.beta       = 2.0; % 1.0;
% %             s.eta        = 0.2;
%             obj.filter = Filter.create(s);
% % 
%             s.filterType = 'FilterAdjointAndProject';   
%             s.mesh       = obj.mesh;
%             s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.filterStep = 'PDE';
%             s.beta       = 2.0; % 1.0;
% %             s.eta        = 0.2;
%             obj.filterAdjointComp = Filter.create(s);
        end

        function createFilterConnectivity(obj)
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filterConnect = f;
% 
%             s.filterType = 'FilterAndProject';
% %             s.filterType = 'CloseOperator'; %'FilterAndProject';
%             s.mesh       = obj.mesh;
%             s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.filterStep = 'PDE';
%             s.beta       = 4.0; % 1.0;
% %             s.eta        = 0.2;
%             obj.filterConnect = Filter.create(s);
% % 
%             s.filterType = 'FilterAdjointAndProject';   
% %             s.filterType = 'CloseAdjointOperator'; %'FilterAdjointAndProject';   
%             s.mesh       = obj.mesh;
%             s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.filterStep = 'PDE';
%             s.beta       = 4.0; % 1.0;
% %             s.eta        = 0.2;
%             obj.filterAdjointConnect = Filter.create(s);
        end

        function createMaterialInterpolator(obj)
            E0   = 1e-5;
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
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filter;
            s.filterAdjoint              = obj.filterAdjointComp;
            s.complainceFromConstitutive = obj.createComplianceFromConstitutive();
            s.material                   = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;
%             v = VolumeConstraint(s);
            s.filter       = obj.filter;
            s.filterAdjoint = obj.filterAdjointComp;
            v = VolumeConstraintWithProjection(s);
            obj.volume = v;
        end

        function createEigenValueConstraint(obj)                           
            s.mesh              = obj.mesh;
            s.designVariable    = obj.designVariable;
            s.filter            = obj.filterConnect;
            s.filterAdjoint     = obj.filterAdjointConnect;
            s.boundaryConditions = obj.createEigenvalueBoundaryConditions();
            s.targetEigenValue  = obj.lambda1min;    
            s.isCompl           = true;
            obj.minimumEigenValue = StiffnessEigenModesConstraint(s);
        end

        function createEigenValue(obj)                           
            s.mesh              = obj.mesh;
            s.designVariable    = obj.designVariable;
            s.filter            = obj.filter;     
            s.boundaryConditions= obj.createEigenvalueBoundaryConditions();
            s.isCompl           = true;
            obj.eigenvalue = MaximumEigenValueFunctional(s);
        end

        function createPerimeter(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filter;
            p = SimplePerimeterFunctional(s);
            obj.perimeter = p;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.shapeFunctions{2} = obj.perimeter;
%             s.shapeFunctions{3} = obj.eigenvalue;
            s.weights           = [1.0,2.5]; %0.5,0.0,1.0v,1.0
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
            s.shapeFunctions{2} = obj.minimumEigenValue; 
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createPrimalUpdater(obj)
            s.mesh = obj.mesh;
            obj.primalUpdater = SLERP(s);
        end

        function createOptimizer(obj)
            s.monitoring       = true;
            s.cost             = obj.cost;
            s.constraint       = obj.constraint;
            s.designVariable   = obj.designVariable;
%             s.GIFname        = '1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+string(obj.eta)+'GIF';
            s.GIFname           = '1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+'GIF';
            s.maxIter           = 500;
            s.tolerance         = 1e-3;
            s.constraintCase{1} = 'EQUALITY';
            s.constraintCase{2} = 'INEQUALITY';      
            s.primalUpdater     = obj.primalUpdater;
            s.etaNorm           = 0.02; 
            s.etaNormMin        = 0.02;
            s.gJFlowRatio       = 0.2; %obj.gJ;
            s.etaMax            = 0.6; %0.1;   %1.0 0.2
            s.etaMaxMin         = 0.02; %0.05; %0.01;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
%             saveas(figure(1),'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+'design'+string(obj.eta)+'.png','png')
%             saveas(figure(2),'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+'graficos'+string(obj.eta)+'.png','png')
%             writematrix(obj.designVariable.fun.fValues,'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+string(obj.eta)+'.txt')
            saveas(figure(1),'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+'design2.png','png')
            saveas(figure(2),'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+'graficos2.png','png')
            writematrix(obj.designVariable.fun.fValues,'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+'.txt')
        end

        function m = createMaterial(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = obj.filter.compute(f{1},1);            
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            s.mesh                 = obj.mesh;
            m = Material.create(s);
        end

        function bc = createElasticBoundaryConditions(obj)
            type = 'cantilever';%' 'bridge'; %  
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