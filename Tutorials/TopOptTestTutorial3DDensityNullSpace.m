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
            %obj.mesh = HexaMesh(1,1,1,20,20,20);
            

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



            % GiD
            file       = 'HexaFine';
            a.fileName = file;
            s          = FemDataContainer(a);
            obj.mesh   = s.mesh;

            %obj.mesh = HexaMesh(1,1,1,10,10,10);
        end

        function createDesignVariable(obj)
            s.type             = 'Full'; % Custom
            g                  = GeometricalFunction(s);
            lsFun              = g.computeLevelSetFunction(obj.mesh);
            s.mesh             = obj.mesh;
            s.order            = 'P1';
            s.fValues          = 1-heaviside(lsFun.fValues);
            s.fun              = LagrangianFunction(s);
            s.type             = 'Density';
            s.plotting         = false;
            dens               = DesignVariable.create(s);
            obj.designVariable = dens;
        end

        function createFilter(obj)
%             s.filterType = 'P1';
%             s.mesh  = obj.mesh;
%             s.test  = LagrangianFunction.create(obj.mesh,1,'P0');
%             s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
%             f = Filter.create(s);
%             obj.filter = f;

            s.filterType   = 'PDE';
            s.mesh         = obj.mesh;
            s.boundaryType = 'Neumann';
            s.metric       = 'Isotropy';
            s.trial        = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.filter     = Filter.create(s);
            epsilon        = 4*obj.mesh.computeMeanCellSize(); % aquest 1 potser el toquem; és el radi del filtre que penalitza tant a l'interior com a la boundary
            obj.filter.updateEpsilon(epsilon);
        end

        function createMaterialInterpolator(obj)
            E0 = 1e-3;
            nu0 = 1/3;
            ndim = obj.mesh.ndim;
            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);


            E1 = 1;
            nu1 = 0.45;
            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.interpolation  = 'SIMP_P3'; % SIMPALL, SIMP_P3
            s.dim            = '3D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;
        end

        function m = createMaterial(obj)
            f = obj.designVariable.fun;
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
            s.volumeTarget = 0.6;
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
            n = obj.mesh.nnodes;
            h = obj.mesh.computeMinCellSize();
            M = h^2*sparse(1:n,1:n,ones(1,n),n,n);
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
            s.maxIter        = 1000;%1000
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.primal         = 'PROJECTED GRADIENT';
            s.ub             = 1;
            s.lb             = 0;
            s.etaNorm        = 0.02;
            s.gJFlowRatio    = 0.7;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
            %obj.designVariable.fun.print('densityNullSpaceLastIter','Paraview');
        end

        function bc = createBoundaryConditionsCube(obj)
            xMax = max(obj.mesh.coord(:,1));
            yMax = max(obj.mesh.coord(:,2));
            zMax = max(obj.mesh.coord(:,3));
            loadXY  = @(coor) coor(:,1)>=0.1*xMax & coor(:,1)<=0.9*xMax & coor(:,2)>=0.1*yMax & coor(:,2)<=0.9*yMax;
            loadYZ  = @(coor) coor(:,3)>=0.1*zMax & coor(:,3)<=0.9*zMax & coor(:,2)>=0.1*yMax & coor(:,2)<=0.9*yMax;
            loadXZ  = @(coor) coor(:,1)>=0.1*xMax & coor(:,1)<=0.9*xMax & coor(:,3)>=0.1*zMax & coor(:,3)<=0.9*zMax;
            isDir   = @(coor)  abs(coor(:,1))==0; % 4 potes
            isForceXu = @(coor)  abs(coor(:,1))==xMax;
            isForceYu = @(coor)  abs(coor(:,2))==yMax & loadXZ(coor);
            isForceYd = @(coor)  abs(coor(:,2))==0 & loadXZ(coor);
            isForceZu = @(coor)  abs(coor(:,3))==zMax & loadXY(coor);
            isForceZd = @(coor)  abs(coor(:,3))==0 & loadXY(coor);

            sDir{1}.domain    = @(coor) isForceYu(coor) | isForceYd(coor) | isForceZu(coor) | isForceZd(coor);
            sDir{1}.direction = [1,2,3];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForceXu(coor) & loadYZ(coor);
            sPL{1}.direction = 1;
            sPL{1}.value     = -1;

            sPL{2}.domain    = @(coor) isDir(coor) & loadYZ(coor);
            sPL{2}.direction = 1;
            sPL{2}.value     = 1;
% 
%             sPL{2}.domain    = @(coor) isForceYu(coor);
%             sPL{2}.direction = 2;
%             sPL{2}.value     = -1;
% 
%             sPL{3}.domain    = @(coor) isForceYd(coor);
%             sPL{3}.direction = 2;
%             sPL{3}.value     = 1;
% 
%             sPL{4}.domain    = @(coor) isForceZu(coor);
%             sPL{4}.direction = 3;
%             sPL{4}.value     = -1;
% 
%             sPL{5}.domain    = @(coor) isForceZd(coor);
%             sPL{5}.direction = 3;
%             sPL{5}.value     = 1;



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