classdef TopOptTestTutorial3DLevelSet < handle

    properties (Access = private)
        mesh
        filter
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        volume
        cost
        constraint
        dualVariable
        optimizer
    end

    methods (Access = public)

        function obj = TopOptTestTutorial3DLevelSet()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createDualVariable();
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            obj.mesh = HexaMesh(6,1,1,120,20,20);
        end

        function createDesignVariable(obj)
            % Density:
            %s.fHandle = @(x) ones(size(x(1,:,:)));
            %s.ndimf   = 1;
            %s.mesh    = obj.mesh;
            %aFun      = AnalyticalFunction(s);
            %s.fun     = aFun.project('P1');
            %s.mesh    = obj.mesh;
            %s.type = 'Density';
            %s.plotting = false;
            %dens    = DesignVariable.create(s);
            %obj.designVariable = dens;

            % LevelSet:
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

        function createMaterialInterpolator(obj)
            %DENSITY
            %E0 = 1e-3;
            %nu0 = 1/3;
            %ndim = obj.mesh.ndim;
            %matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            %matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);


            %E1 = 1;
            %nu1 = 1/3;
            %matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            %matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            %s.interpolation  = 'SIMPALL';
            %s.dim            = '3D';
            %s.matA = matA;
            %s.matB = matB;

            %m = MaterialInterpolator.create(s);
            %obj.materialInterpolator = m;

            %LEVEL-SET
            E0   = 1e-3;
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
            s.dim            = '3D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;

        end
  
        function m = createMaterial(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = obj.filter.compute(f,1);            
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '3D';
            m = Material.create(s);
        end
    

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '3D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'rMINRES';
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function c = createComplianceFromConstiutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstiutiveTensor(s);
        end

        function createCompliance(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filter;
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;    %volume target
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
           s.test  = LagrangianFunction.create(obj.mesh,1,'P1');
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh  = obj.mesh;
            s.type  = 'MassMatrix';
            LHS = LHSintegrator.create(s);
            M = LHS.compute;

            h = obj.mesh.computeMinCellSize();
            M = h^2*eye(size(M));
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
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
            s.maxIter        = 250;
            s.tolerance      = 1e-8; 
            s.constraintCase = {'EQUALITY'};
            s.primal         = 'SLERP'; 
            s.ub             = inf;
            s.lb             = -inf;
            s.etaNorm        = 0.02;
            s.gJFlowRatio    = 2; 
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
            obj.designVariable.fun.print('LevelSet_3D_2'); %Guarda la simulació automàticament per poder veure-la després a paraview
        end

        function bc = createBoundaryConditions(obj)
           
%... PER A 2D
           % xMax    = max(obj.mesh.coord(:,1));
            %yMax    = max(obj.mesh.coord(:,2));

            %isDir1   = @(coor)  (abs(coor(:,2))==0 &  abs(coor(:,1))>=0 & abs(coor(:,1))<=0.05*xMax); % is Dirichlet1 = on els desplaçaments estan imposats (banda esquerra)
            %isDir2   = @(coor)  (abs(coor(:,2))==0 &
            %abs(coor(:,1))>=0.95*xMax & abs(coor(:,1))<=xMax); % is Dirichlet2 = on els desplaçaments estan imposats (banda dreta)%

            %isForce1 = @(coor)  (abs(coor(:,2))==yMax & abs(coor(:,1))>=0.475*xMax & abs(coor(:,1))<=0.525*xMax); % isForce1 = força; amunt, centrat


            %sDir{1}.domain    = @(coor) isDir1(coor);
            %sDir{1}.direction = 2;
            %sDir{1}.value     = 0;

            %sDir{2}.domain    = @(coor) isDir2(coor);
            %sDir{2}.direction = [1,2];
            %sDir{2}.value     = 0;

            %sPL{1}.domain    = @(coor) isForce1(coor);
            %sPL{1}.direction = 2;
            %sPL{1}.value     = -1;

%...  PER A 3D  
            xMax = max(obj.mesh.coord(:,1));
            yMax = max(obj.mesh.coord(:,2));
            zMax = max(obj.mesh.coord(:,3));
        
            isDir1   = @(coor)  (abs(coor(:,2))==0 &  abs(coor(:,1))>=0 & abs(coor(:,1))<=0.05*xMax); % is Dirichlet1 = on els desplaçaments estan imposats (banda esquerra)
            isDir2   = @(coor)  (abs(coor(:,2))==0 & abs(coor(:,1))>=0.95*xMax & abs(coor(:,1))<=xMax); % is Dirichlet2 = on els desplaçaments estan imposats (banda dreta)

            isForce1 = @(coor)  (abs(coor(:,2))==yMax & abs(coor(:,1))>=0.475*xMax & abs(coor(:,1))<=0.525*xMax & abs(coor(:,3))>=0.4*zMax & abs(coor(:,3))<=0.6*zMax); % isForce1 = força; amunt, centrat

            sDir{1}.domain    = @(coor) isDir1(coor);
            sDir{1}.direction = [2,3];
            sDir{1}.value     = 0;
            
            sDir{2}.domain    = @(coor) isDir2(coor);
            sDir{2}.direction = [1,2,3];
            sDir{2}.value     = 0;

            sPL{1}.domain    = @(coor) isForce1(coor);
            sPL{1}.direction = 3;
            sPL{1}.value     = -1;

            %... AIXÒ ES QUEDA IGUAL
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