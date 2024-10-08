classdef BulkCaseStudy < handle

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

        function obj = BulkCaseStudy()
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
            obj.mesh = TetraMesh(1,1,1,20,20,20);
        end

        function createDesignVariable(obj)
            s.type             = 'Custom';
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
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
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
            s.boundaryConditions = obj.createBoundaryConditions();
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
            s.maxIter        = 10000;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.ub             = 1;
            s.lb             = 0;
            s.primal         = 'PROJECTED GRADIENT';
            opt = OptimizerMMA(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj)
            xMax = max(obj.mesh.coord(:,1));
            yMax = max(obj.mesh.coord(:,2));
            zMax = max(obj.mesh.coord(:,3));

            % Edges of face z=0
            edge1 = @(coor) coor(:,3)==0 & coor(:,2)==0;
            edge2 = @(coor) coor(:,3)==0 & coor(:,1)==xMax;
            edge3 = @(coor) coor(:,3)==0 & coor(:,2)==yMax;
            edge4 = @(coor) coor(:,3)==0 & coor(:,1)==0;

            % Edges of face z=zMax
            edge5 = @(coor) coor(:,3)==zMax & coor(:,2)==0;
            edge6 = @(coor) coor(:,3)==zMax & coor(:,1)==xMax;
            edge7 = @(coor) coor(:,3)==zMax & coor(:,2)==yMax;
            edge8 = @(coor) coor(:,3)==zMax & coor(:,1)==0;

            % Edges connecting faces z=0 with z=zMax
            edge9  = @(coor) coor(:,2)==0 & coor(:,1)==0;
            edge10 = @(coor) coor(:,2)==0 & coor(:,1)==xMax;
            edge11 = @(coor) coor(:,2)==yMax & coor(:,1)==xMax;
            edge12 = @(coor) coor(:,2)==yMax & coor(:,1)==0;

            isDir   = @(coor)  edge1(coor) | edge2(coor) | edge3(coor) | edge4(coor) ...
                | edge5(coor) | edge6(coor) | edge7(coor) | edge8(coor) | edge9(coor)...
                | edge10(coor) | edge11(coor) | edge12(coor);

            isForceXu = @(coor) coor(:,1)==xMax;
            isForceXd = @(coor) coor(:,1)==0;
            isForceYu = @(coor) coor(:,2)==yMax;
            isForceYd = @(coor) coor(:,2)==0;
            isForceZu = @(coor) coor(:,3)==zMax;
            isForceZd = @(coor) coor(:,3)==0;

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2,3];
            sDir{1}.value     = 0;

            delta = 0.15;

            sPL{1}.domain    = @(coor) isForceXu(coor);
            sPL{1}.direction = 1;
            sPL{1}.value     = 1;

            sPL{2}.domain    = @(coor) isForceXd(coor);
            sPL{2}.direction = 1;
            sPL{2}.value     = -(1+delta);

            sPL{3}.domain    = @(coor) isForceYu(coor);
            sPL{3}.direction = 2;
            sPL{3}.value     = 1;

            sPL{4}.domain    = @(coor) isForceYd(coor);
            sPL{4}.direction = 2;
            sPL{4}.value     = -(1+delta);

            sPL{5}.domain    = @(coor) isForceZu(coor);
            sPL{5}.direction = 3;
            sPL{5}.value     = 1;

            sPL{6}.domain    = @(coor) isForceZd(coor);
            sPL{6}.direction = 3;
            sPL{6}.value     = -(1+delta);

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