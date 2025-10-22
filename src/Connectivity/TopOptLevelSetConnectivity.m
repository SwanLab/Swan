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
        filterAdjointProj 
        filterProj
        complianceProj
        volumeProj
        eta
        beta
        primalUpdater
        eigenvalue
        c
        type
        p
    end 

    methods (Access = public)
        function obj = TopOptLevelSetConnectivity()
            for type = ["cantilever"]
%                 for p = [0.0]
                for c = [2,1,4] %0.15]0.05,0.15,0.5,1,3,4,"cantilever",
                    for lambda1min = [0.6] %0.15]0.05,0.15,0.5,
                        obj.c = c;
                        obj.p = 0.0;
                        obj.type = type;
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
                        obj.createConductivityInterpolator();
                        obj.createMassInterpolator();
                        obj.createEigenValueConstraint();   
                        obj.createEigenValue()          
                        obj.createPerimeter();                  
                        obj.createVolumeConstraint();
                        obj.createCost();
                        obj.createConstraint();
                        obj.createPrimalUpdater();
                        obj.createOptimizer();
                    end
                end
            end
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
            x1      = linspace(0,2,100);
            x2      = linspace(0,1,50);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh.create(s);
        end

        function createDesignVariable(obj)
            s.type = 'Full';
            g      = GeometricalFunction(s);
            lsFun  = g.computeLevelSetFunction(obj.mesh);
%             lsFun = LagrangianFunction.create(obj.mesh,1,'P1');
%             lsFun.setFValues(importdata('_cantilever_case_5_p_00.8.txt'))
            s.fun  = lsFun;
            s.mesh = obj.mesh;
            s.type = 'LevelSet';
            s.plotting = true;
            ls     = DesignVariable.create(s);
            obj.designVariable = ls;
            ls.fun.print('teste2D','Paraview')
        end

        function createFilter(obj)
            s.filterType = 'PDE';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;

            if isequal(obj.c, 5)
                s.filterType = 'PDE';
                s.mesh       = obj.mesh;
                s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
                f            = Filter.create(s);
                obj.filter = f;
            end
        end

        function createFilterConnectivity(obj)
           if isequal(obj.c, 1)
                s.filterType = 'FilterAndProject';
                s.mesh       = obj.mesh;
                s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
                s.filterStep = 'PDE';
                s.beta       = 20.0; 
                obj.filterConnect = Filter.create(s);
     
                s.filterType = 'FilterAdjointAndProject';   
                s.mesh       = obj.mesh;
                s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
                s.filterStep = 'PDE';
                s.beta       = 20.0; 
                obj.filterAdjointConnect = Filter.create(s);
           elseif isequal(obj.c, 2)
                s.filterType = 'FilterAndProject';
                s.mesh       = obj.mesh;
                s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
                s.filterStep = 'PDE';
                s.beta       = 4.0; 
                s.eta        = 0.2;
                obj.filterConnect = Filter.create(s);

                s.filterType = 'FilterAdjointAndProject';   
                s.mesh       = obj.mesh;
                s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
                s.filterStep = 'PDE';
                s.beta       = 4.0; 
                s.eta        = 0.2;
                obj.filterAdjointConnect = Filter.create(s);

           elseif isequal(obj.c, 3)
                s.filterType = 'PDE';
                s.mesh       = obj.mesh;
                s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
                f            = Filter.create(s);
                obj.filterConnect = f;
                obj.filterAdjointConnect =[];

          elseif isequal(obj.c, 4)
                s.filterType = 'FilterAndProject';
                s.mesh       = obj.mesh;
                s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
                s.filterStep = 'PDE';
                s.beta       = 4.0; % 20.0
                s.eta        = 0.5;
                obj.filterConnect = Filter.create(s);
     
                s.filterType = 'FilterAdjointAndProject';   
                s.mesh       = obj.mesh;
                s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
                s.filterStep = 'PDE';
                s.beta       = 4.0;  %20.0
                s.eta        = 0.5;
                obj.filterAdjointConnect = Filter.create(s);
          elseif isequal(obj.c, 5)
                s.filterType = 'LUMP';
                s.mesh       = obj.mesh;
                s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
                f            = Filter.create(s);
                obj.filterConnect = f;
                obj.filterAdjointConnect =[];
            end
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
            s.filter = obj.filter;
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

        function createEigenValue(obj)                           
            s.mesh              = obj.mesh;
            s.designVariable    = obj.designVariable;
            s.filter            = obj.filterConnect; 
            s.filterAdjoint     = obj.filterAdjointConnect;
            s.conductivityInterpolator = obj.conductivityInterpolator; 
            s.massInterpolator         = obj.massInterpolator; 
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
%             s.shapeFunctions{2} = obj.eigenvalue;
            s.weights           = [1.0,obj.p]; %0.5,0.0,1.0v,1.0,5.0,3.0 ,2.5
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
            s.GIFname           = 'beta20_'+string(obj.type)+'_case_'+string(obj.c)+'_p_'+string(obj.p)+'_lamb_'+string(obj.lambda1min)+'-';
            s.maxIter           = 1500;
            s.tolerance         = 1e-3;
            s.constraintCase{1} = 'EQUALITY';
            s.constraintCase{2} = 'INEQUALITY';      
            s.primalUpdater     = obj.primalUpdater;
            s.etaNorm           = 0.02; 
            s.etaNormMin        = 0.02;
%             s.gJFlowRatio       = 2.0; %obj.gJ;
%             s.etaMax            = 0.6; %0.1;   %1.0 0.2
%             s.etaMaxMin         = 0.02; %0.05; %0.01;
            s.gJFlowRatio       = 2.0; %obj.gJ;
            s.etaMax            = 1.0; %1.0; %0.1;   %1.0 0.2
            s.etaMaxMin         = 0.02; %0.05; %0.01;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
%             saveas(figure(1),'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+'design'+string(obj.eta)+'.png','png')
%             saveas(figure(2),'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+'graficos'+string(obj.eta)+'.png','png')
%             writematrix(obj.designVariable.fun.fValues,'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+string(obj.eta)+'.txt')
            saveas(figure(1),'beta20_'+'_'+string(obj.type)+'_case_'+string(obj.c)+'_p_'+string(obj.p)+string(obj.lambda1min)+'design2.png','png')
            saveas(figure(2),'beta20_'+'_'+string(obj.type)+'_case_'+string(obj.c)+'_p_'+string(obj.p)+string(obj.lambda1min)+'graficos2.png','png')
            writematrix(obj.designVariable.fun.fValues,'beta20_'+'_'+string(obj.type)+'_case_'+string(obj.c)+'_p_'+string(obj.p)+string(obj.lambda1min)+'.txt')
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
%             type = 'cantilever';%' 'bridge'; %  
            type = obj.type;
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