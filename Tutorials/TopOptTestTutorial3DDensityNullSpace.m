classdef TopOptTestTutorial3DDensityNullSpace < handle

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

        function obj = TopOptTestTutorial3DDensityNullSpace()
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

            %Per cub:
            obj.mesh = HexaMesh(1,1,1,20,20,20);

            % Cas cilindre linux:
%             m2D = QuadMesh(1,1,20,20);
%             s.type        = 'Circle';
%             s.radius      = 0.15;
%             s.xCoorCenter = 0.5;
%             s.yCoorCenter = 0.5;
%             gFun          = GeometricalFunction(s);
%             ls            = gFun.computeLevelSetFunction(m2D);
%             sUm.backgroundMesh = m2D;
%             sUm.boundaryMesh   = m2D.createBoundaryMesh();
%             uMesh = UnfittedMesh(sUm);
%             uMesh.compute(ls.fValues);
%             IM = uMesh.createInnerMesh();
%             meshCylinder = IM.provideExtrudedMesh(1);

            % Cas cilindre windows:
%             load('meshCylinder.mat','meshCylinder');
%             obj.mesh = meshCylinder;
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);
            s.fun     = aFun.project('P1');
            s.mesh    = obj.mesh;
            s.type = 'Density';
            s.plotting = false;
            dens    = DesignVariable.create(s);
            obj.designVariable = dens;
        end

        function createFilter(obj)
            s.filterType = 'P1';
            s.mesh  = obj.mesh;
            s.test  = LagrangianFunction.create(obj.mesh,1,'P0');
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;

%             s.filterType   = 'PDE';
%             s.mesh         = obj.mesh;
%             s.boundaryType = 'Robin';
%             s.metric       = 'Isotropy';
%             s.trial        = LagrangianFunction.create(obj.mesh,1,'P1');
%             obj.filter     = Filter.create(s);
%             epsilon        = 1*obj.mesh.computeMeanCellSize(); % aquest 1 potser el toquem; és el radi del filtre que penalitza tant a l'interior com a la boundary
%             obj.filter.updateEpsilon(epsilon);
        end

        function createMaterialInterpolator(obj)
            E0 = 1e-3;
            nu0 = 1/3;
            ndim = obj.mesh.ndim;
            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);


            E1 = 1;
            nu1 = 0.499;
            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

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
            f = f.project('P1');            
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
            s.boundaryConditions = obj.createBoundaryConditionsCube(); % obj.createBoundaryConditionsCube      //  obj.createBoundaryConditionsCylinder
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
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
            s.volumeTarget = 0.7;
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
            s.maxIter        = 700;%1000
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.primal         = 'PROJECTED GRADIENT';
            s.ub             = 1;
            s.lb             = 0;
            s.etaNorm        = 0.05;
            s.gJFlowRatio    = 0.6;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditionsCube(obj)
            xMax = max(obj.mesh.coord(:,1));
            yMax = max(obj.mesh.coord(:,2));
            zMax = max(obj.mesh.coord(:,3));
            isDir   = @(coor)  abs(coor(:,1))==0; % 4 potes
            isForceXu = @(coor)  abs(coor(:,1))==xMax;
            isForceYu = @(coor)  abs(coor(:,2))==yMax;
            isForceYd = @(coor)  abs(coor(:,2))==0;
            isForceZu = @(coor)  abs(coor(:,3))==zMax;
            isForceZd = @(coor)  abs(coor(:,3))==0;

            sDir{1}.domain    = @(coor) isDir(coor) | isForceYu(coor) | isForceYd(coor) | isForceZu(coor) | isForceZd(coor);
            sDir{1}.direction = [1,2,3];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForceXu(coor);
            sPL{1}.direction = 1; % direccio x +
            sPL{1}.value     = -1; % sentit -

  

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


        function bc = createBoundaryConditionsCylinder(obj)
            zMax = max(obj.mesh.coord(:,3));
            isDir     = @(coor)  sqrt((abs(coor(:,1))-0.5).^2+(abs(coor(:,2))-0.5).^2)>=0.148; % 0.148 ??
            isForceZu = @(coor)  abs(coor(:,3))==zMax;
            isForceZd = @(coor)  abs(coor(:,3))==0;
% 
            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2,3];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForceZu(coor);
            sPL{1}.direction = 3;
            sPL{1}.value     = -1;

            sPL{2}.domain    = @(coor) isForceZd(coor);
            sPL{2}.direction = 3;
            sPL{2}.value     = 1;

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