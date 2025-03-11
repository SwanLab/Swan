classdef TopOptLevelSet3DConnectivity< handle

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
        gJ
        lambda1min
        filterAdjointProj 
        filterProj
        complianceProj
        volumeProj
        eta
    end 

    methods (Access = public)
        function obj = TopOptLevelSet3DConnectivity()
            for gJ = [2.0]
                for lambda1min = [0.05]
%                 for eta = [1.0]
                    obj.gJ = gJ;
                    obj.lambda1min = lambda1min;
%                     obj.eta        = eta;
                    obj.init()
                    obj.createMesh();
                    obj.createDesignVariable();
                    obj.createFilter();
                    obj.createFilterCompliance();
                    obj.createFilterComplianceProjected();
                    obj.createFilterConnectivity();
                    obj.createMaterialInterpolator();
                    obj.createElasticProblem();
                    obj.createComplianceFromConstitutive();
                    obj.createCompliance();
                    obj.createComplianceProjected();
                    obj.createEigenValueConstraint();                             
                    obj.createVolumeConstraint();
                    obj.createVolumeConstraintProjected();
                    obj.createCost();
                    obj.createConstraint();
                    obj.createDualVariable();
                    obj.createOptimizer();
%                 end
                end
            end
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            %UnitMesh better
            obj.mesh = HexaMesh(4,4,3,20,20,15);
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

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

         function createFilterCompliance(obj)
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filterComp = f;
         end

         function createFilterComplianceProjected(obj)
            s.filterType = 'FilterAndProject';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            s.filterStep = 'LUMP';
            s.beta       = 2.0;
            s.eta        = obj.eta;
            obj.filterProj = Filter.create(s);

            s.filterType = 'FilterAdjointAndProject';   
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            s.filterStep = 'LUMP';
            s.beta       = 2.0;
            s.eta        = obj.eta;
            obj.filterAdjointProj = Filter.create(s);
         end 

        function createComplianceProjected(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filterProj;
            if ~isempty(obj.filterAdjointProj)
                s.filterAdjoint               = obj.filterAdjointProj;
            end
            s.complainceFromConstitutive  = obj.createComplianceFromConstitutive();
            s.material                    = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.complianceProj = c;
        end

        function createFilterConnectivity(obj)
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filterConnect = f;

%             s.filterType = 'FilterAndProject';
%             s.mesh       = obj.mesh;
%             s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.filterStep = 'LUMP';
%             s.beta       = 2.0;
%             s.eta        = obj.eta;
%             obj.filterConnect = Filter.create(s);
% 
%             s.filterType = 'FilterAdjointAndProject';   
%             s.mesh       = obj.mesh;
%             s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.filterStep = 'LUMP';
%             s.beta       = 2.0;
%             s.eta        = obj.eta;
%             obj.filterAdjointConnect = Filter.create(s);
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
% 
%             s.typeOfMaterial = 'ISOTROPIC';
%             s.interpolation  = 'SIMPALL';
%             s.dim            = '2D';
%             s.matA = matA;
%             s.matB = matB;
% 
            s.interpolation  = 'SIMP_P3';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '3D';
            s.boundaryConditions = obj.createElasticBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'CG';
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
            s.material                   = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end
       
     
        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;
            v = VolumeConstraint(s);
            obj.volume = v;
        end

       function createVolumeConstraintProjected(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filterProj;
            s.filterAdjoint = obj.filterAdjointProj;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;
            v = VolumeConstraintWithProjection(s);
            obj.volumeProj = v;
        end

        function createEigenValueConstraint(obj)                           
            s.mesh              = obj.mesh;
            s.designVariable    = obj.designVariable;
            s.filter            = obj.filterConnect;
            s.filterAdjoint     = obj.filterAdjointConnect;  
            s.boundaryConditions = obj.createElasticBoundaryConditions();
            s.targetEigenValue  = obj.lambda1min;       
            obj.minimumEigenValue = StiffnessEigenModesConstraint(s);
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
%             s.shapeFunctions{2} = obj.complianceProj;
            s.weights           = [1.0];
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

%         function M = createMassMatrix(obj)
%             s.test  = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.mesh  = obj.mesh;
%             s.type  = 'MassMatrix';
%             LHS = LHSIntegrator.create(s);
%             M = LHS.compute;
% 
%             h = obj.mesh.computeMinCellSize();
%             M = h^2*eye(size(M));
%         end

        function M = createMassMatrix(obj)
            n = obj.mesh.nnodes;
            h = obj.mesh.computeMinCellSize();
            M = h^2*sparse(1:n,1:n,ones(1,n),n,n);
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
%             s.shapeFunctions{2} = obj.minimumEigenValue; 
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createDualVariable(obj)
            s.nConstraints   = 1;
            l                = DualVariable(s);
            obj.dualVariable = l;
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
%             s.GIFname        = '1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+string(obj.eta)+'GIF';
            s.GIFname        = '1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+'GIF';
            s.maxIter        = 1000;
            s.tolerance      = 1e-8;
            s.constraintCase{1} = 'EQUALITY';
%             s.constraintCase{2} = 'INEQUALITY';      
            s.primal         = 'SLERP';
            s.ub             = inf;
            s.lb             = -inf;
            s.etaNorm        = 0.03; %0.02 0.1
            s.etaNormMin     = 0.02;
            s.gJFlowRatio    = obj.gJ;
            s.etaMax         = 1.0;   %1.0
            s.etaMaxMin      = 0.01; %0.01;
            s.filter         = obj.filterComp;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
%             saveas(figure(1),'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+'design'+string(obj.eta)+'.png','png')
%             saveas(figure(2),'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+'graficos'+string(obj.eta)+'.png','png')
%             writematrix(obj.designVariable.fun.fValues,'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+string(obj.eta)+'.txt')
            saveas(figure(1),'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+'design.png','png')
            saveas(figure(2),'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+'graficos.png','png')
            writematrix(obj.designVariable.fun.fValues,'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+'.txt')
        end

        function m = createMaterial(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = obj.filterComp.compute(f{1},1);            
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '3D';
            s.mesh                 = obj.mesh;
            m = Material.create(s);
        end

        function bc = createElasticBoundaryConditions(obj)
            xMax = max(obj.mesh.coord(:,1));
            yMax = max(obj.mesh.coord(:,2));
            zMax = max(obj.mesh.coord(:,3));
            isDir   = @(coor)  (abs(coor(:,3))==0 & abs(coor(:,2))>=0.4*yMax & abs(coor(:,2))<=0.6*yMax & abs(coor(:,1))>=0.4*xMax & abs(coor(:,1))<=0.6*xMax);
            isForce = @(coor)  (abs(coor(:,3))==zMax);

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2,3];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 3;
            sPL{1}.value     = -1;

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
            zMin    = min(obj.mesh.coord(:,3));
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            zMax    = max(obj.mesh.coord(:,3));
            isDown  = @(coor) abs(coor(:,2))==zMin;
            isUp    = @(coor) abs(coor(:,2))==zMax;
            isLeft  = @(coor) abs(coor(:,1))==xMin;
            isRight = @(coor) abs(coor(:,1))==xMax;
            isFront = @(coor) abs(coor(:,1))==yMin;
            isBack = @(coor) abs(coor(:,1))== yMax;
            isDir   = @(coor)  isDown(coor) | isUp(coor) | isLeft(coor) | isRight(coor) | isFront(coor) | isBack(coor);  
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