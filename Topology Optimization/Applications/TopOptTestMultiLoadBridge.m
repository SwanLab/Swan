classdef TopOptTestMultiLoadBridge < handle

    properties (Access = private)
        mesh
        filter
        designVariable
        materialInterpolator
        physicalProblemLeft
        physicalProblemCenter
        physicalProblemRight
        complianceLeft
        complianceCenter
        complianceRight
        volume
        cost
        constraint
        dualVariable
        optimizer
    end

    methods (Access = public)

        function obj = TopOptTestMultiLoadBridge()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblemLeft();
            obj.createElasticProblemCenter();
            obj.createElasticProblemRight();
            obj.createComplianceLeft();
            obj.createComplianceCenter();
            obj.createComplianceRight();
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
            x1       = linspace(0,10,200);
            x2       = linspace(0,2,40);
            [xv,yv]  = meshgrid(x1,x2);
            [F,V]    = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh.create(s);
        end

        function createDesignVariable(obj)
            s.type             = 'Full';
            g                  = GeometricalFunction(s);
            lsFun              = g.computeLevelSetFunction(obj.mesh);
            s.fun              = lsFun;
            s.mesh             = obj.mesh;
            s.type             = 'LevelSet';
            s.plotting         = true;
            ls                 = DesignVariable.create(s);
            obj.designVariable = ls;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filter   = f;
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

            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;
        end

        function createElasticProblemLeft(obj)
            s.mesh                  = obj.mesh;
            s.scale                 = 'MACRO';
            s.material              = obj.createMaterial();
            s.dim                   = '2D';
            s.boundaryConditions    = obj.createBoundaryConditions(0,1/9);
            s.interpolationType     = 'LINEAR';
            s.solverType            = 'REDUCED';
            s.solverMode            = 'DISP';
            fem                     = ElasticProblem(s);
            obj.physicalProblemLeft = fem;
        end

        function createElasticProblemCenter(obj)
            s.mesh                    = obj.mesh;
            s.scale                   = 'MACRO';
            s.material                = obj.createMaterial();
            s.dim                     = '2D';
            s.boundaryConditions      = obj.createBoundaryConditions(4/9,5/9);
            s.interpolationType       = 'LINEAR';
            s.solverType              = 'REDUCED';
            s.solverMode              = 'DISP';
            fem                       = ElasticProblem(s);
            obj.physicalProblemCenter = fem;
        end

        function createElasticProblemRight(obj)
            s.mesh                   = obj.mesh;
            s.scale                  = 'MACRO';
            s.material               = obj.createMaterial();
            s.dim                    = '2D';
            s.boundaryConditions     = obj.createBoundaryConditions(8/9,1);
            s.interpolationType      = 'LINEAR';
            s.solverType             = 'REDUCED';
            s.solverMode             = 'DISP';
            fem                      = ElasticProblem(s);
            obj.physicalProblemRight = fem;
        end

        function c = createComplianceFromConstiutive(obj,physicalProblem)
            s.mesh         = obj.mesh;
            s.stateProblem = physicalProblem;
            c              = ComplianceFromConstiutiveTensor(s);
        end

        function createComplianceLeft(obj)
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filter;
            s.complainceFromConstitutive = obj.createComplianceFromConstiutive(obj.physicalProblemLeft);
            s.material                   = obj.createMaterial();
            c                            = ComplianceFunctional(s);
            obj.complianceLeft           = c;
        end

        function createComplianceCenter(obj)
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filter;
            s.complainceFromConstitutive = obj.createComplianceFromConstiutive(obj.physicalProblemCenter);
            s.material                   = obj.createMaterial();
            c                            = ComplianceFunctional(s);
            obj.complianceCenter         = c;
        end

        function createComplianceRight(obj)
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filter;
            s.complainceFromConstitutive = obj.createComplianceFromConstiutive(obj.physicalProblemRight);
            s.material                   = obj.createMaterial();
            c                            = ComplianceFunctional(s);
            obj.complianceRight          = c;
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.complianceLeft;
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
            s.maxIter        = 1000;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.volumeTarget   = 0.4;
            s.primal         = 'SLERP';
            s.etaNorm        = 0.02;
            s.gJFlowRatio    = 1.5;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function m = createMaterial(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = obj.filter.compute(f,'LINEAR');            
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            m = Material.create(s);
        end

        function bc = createBoundaryConditions(obj,t1,t2)
            xMax     = max(obj.mesh.coord(:,1));
            yMax     = max(obj.mesh.coord(:,2));
            isDirX   = @(coor) abs(coor(:,1))<=0.1*xMax | abs(coor(:,1))>=0.9*xMax;
            isDirY   = @(coor) abs(coor(:,2))==0;
            isDir    = @(coor) isDirX(coor) & isDirY(coor);
            isForceX = @(coor) abs(coor(:,1))>=t1*xMax & abs(coor(:,1))<=t2*xMax;
            isForceY = @(coor) abs(coor(:,2))==yMax;
            isForce  = @(coor) isForceX(coor) & isForceY(coor);

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
                pl = PointLoad(obj.mesh, sPL{i});
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);
        end
    end
end