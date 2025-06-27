classdef TopOptTestMultiLoadBridgeAugmented < handle

    properties (Access = private)
        mesh
        filter
        designVariable
        designVariableTar
        materialInterpolator
        physicalProblemLeft
        physicalProblem12
        physicalProblem23
        physicalProblem34
        physicalProblemCenter
        physicalProblem56
        physicalProblem67
        physicalProblem78
        physicalProblemRight
        targetCompliance
        complianceLeft
        compliance12
        compliance23
        compliance34
        complianceCenter
        compliance56
        compliance67
        compliance78
        complianceRight
        volume
        cost
        constraint
        dualVariable
        optimizer
        nLoads
        rho
    end

    methods (Access = public)

        function obj = TopOptTestMultiLoadBridgeAugmented(nL,rho)
            obj.init(nL,rho)
            obj.createMesh();
            obj.createDesignVariable();
            obj.createDesignVariableTargetCompliance();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblemLeft();
            obj.createElasticProblem12();
            obj.createElasticProblem23();
            obj.createElasticProblem34();
            obj.createElasticProblemCenter();
            obj.createElasticProblem56();
            obj.createElasticProblem67();
            obj.createElasticProblem78();
            obj.createElasticProblemRight();
            obj.computeTargetCompliance();
            obj.createComplianceLeft();
            obj.createCompliance12();
            obj.createCompliance23();
            obj.createCompliance34();
            obj.createComplianceCenter();
            obj.createCompliance56();
            obj.createCompliance67();
            obj.createCompliance78();
            obj.createComplianceRight();
            obj.createVolumeFunctional();
            obj.createCost();
            obj.createConstraint();
            obj.createDualVariable();
            obj.createOptimizer();

            saveas(gcf,['NullSLERPResults/TopOpt/MultiLoadBridge/Augmented/Monitoring_rho',num2str(obj.rho),'9Loads.fig']);
            obj.designVariable.fun.print(['NullSLERPResults/TopOpt/MultiLoadBridge/Augmented/rho',num2str(obj.rho),'9Loads_fValues']);
        end

    end

    methods (Access = private)

        function init(obj,nL,r)
            close all;
            obj.nLoads = nL;
            obj.rho = r;
        end

        function createMesh(obj)
            x1       = linspace(0,10,400);
            x2       = linspace(0,2,80);
            [xv,yv]  = meshgrid(x1,x2);
            [F,V]    = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh.create(s);
        end

        function createDesignVariable(obj)
            s.type             = 'Holes';
            s.dim              = 2;
            s.nHoles           = [4,3];
            s.totalLengths     = [2,1];
            s.phiZero          = 0.4;
            s.phases           = [0,0];
            g                  = GeometricalFunction(s);
            lsFun              = g.computeLevelSetFunction(obj.mesh);
            lsFun.fValues      = lsFun.fValues-0.5;
            s.fun              = lsFun;
            s.mesh             = obj.mesh;
            s.type             = 'LevelSet';
            s.plotting         = false;
            ls                 = DesignVariable.create(s);
            obj.designVariable = ls;
        end

        function createDesignVariableTargetCompliance(obj)
            s.type             = 'Holes';
            s.dim              = 2;
            s.nHoles           = [4,3];
            s.totalLengths     = [2,1];
            s.phiZero          = 0.4;
            s.phases           = [0,0];
            g                  = GeometricalFunction(s);
            lsFun              = g.computeLevelSetFunction(obj.mesh);
            s.fun              = lsFun;
            s.mesh             = obj.mesh;
            s.type             = 'LevelSet';
            s.plotting         = false;
            ls                 = DesignVariable.create(s);
            obj.designVariableTar = ls;
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
            x                       = obj.designVariable;
            s.mesh                  = obj.mesh;
            s.scale                 = 'MACRO';
            s.material              = obj.createMaterial(x);
            s.dim                   = '2D';
            s.boundaryConditions    = obj.createBoundaryConditions(0,1/9);
            s.interpolationType     = 'LINEAR';
            s.solverType            = 'REDUCED';
            s.solverMode            = 'DISP';
            s.solverCase = 'rMINRES';
            fem                     = ElasticProblem(s);
            obj.physicalProblemLeft = fem;
        end

        function createElasticProblem12(obj)
            x                       = obj.designVariable;
            s.mesh                  = obj.mesh;
            s.scale                 = 'MACRO';
            s.material              = obj.createMaterial(x);
            s.dim                   = '2D';
            s.boundaryConditions    = obj.createBoundaryConditions(1/9,2/9);
            s.interpolationType     = 'LINEAR';
            s.solverType            = 'REDUCED';
            s.solverMode            = 'DISP';
            s.solverCase = 'rMINRES';
            fem                     = ElasticProblem(s);
            obj.physicalProblem12 = fem;
        end

        function createElasticProblem23(obj)
            x                       = obj.designVariable;
            s.mesh                  = obj.mesh;
            s.scale                 = 'MACRO';
            s.material              = obj.createMaterial(x);
            s.dim                   = '2D';
            s.boundaryConditions    = obj.createBoundaryConditions(2/9,3/9);
            s.interpolationType     = 'LINEAR';
            s.solverType            = 'REDUCED';
            s.solverMode            = 'DISP';
            s.solverCase = 'rMINRES';
            fem                     = ElasticProblem(s);
            obj.physicalProblem23 = fem;
        end

        function createElasticProblem34(obj)
            x                       = obj.designVariable;
            s.mesh                  = obj.mesh;
            s.scale                 = 'MACRO';
            s.material              = obj.createMaterial(x);
            s.dim                   = '2D';
            s.boundaryConditions    = obj.createBoundaryConditions(3/9,4/9);
            s.interpolationType     = 'LINEAR';
            s.solverType            = 'REDUCED';
            s.solverMode            = 'DISP';
            s.solverCase = 'rMINRES';
            fem                     = ElasticProblem(s);
            obj.physicalProblem34 = fem;
        end

        function createElasticProblemCenter(obj)
            x                         = obj.designVariable;
            s.mesh                    = obj.mesh;
            s.scale                   = 'MACRO';
            s.material                = obj.createMaterial(x);
            s.dim                     = '2D';
            s.boundaryConditions      = obj.createBoundaryConditions(4/9,5/9);
            s.interpolationType       = 'LINEAR';
            s.solverType              = 'REDUCED';
            s.solverMode              = 'DISP';
            s.solverCase = 'rMINRES';
            fem                       = ElasticProblem(s);
            obj.physicalProblemCenter = fem;
        end

        function createElasticProblem56(obj)
            x                       = obj.designVariable;
            s.mesh                  = obj.mesh;
            s.scale                 = 'MACRO';
            s.material              = obj.createMaterial(x);
            s.dim                   = '2D';
            s.boundaryConditions    = obj.createBoundaryConditions(5/9,6/9);
            s.interpolationType     = 'LINEAR';
            s.solverType            = 'REDUCED';
            s.solverMode            = 'DISP';
            s.solverCase = 'rMINRES';
            fem                     = ElasticProblem(s);
            obj.physicalProblem56 = fem;
        end

        function createElasticProblem67(obj)
            x                       = obj.designVariable;
            s.mesh                  = obj.mesh;
            s.scale                 = 'MACRO';
            s.material              = obj.createMaterial(x);
            s.dim                   = '2D';
            s.boundaryConditions    = obj.createBoundaryConditions(6/9,7/9);
            s.interpolationType     = 'LINEAR';
            s.solverType            = 'REDUCED';
            s.solverMode            = 'DISP';
            s.solverCase = 'rMINRES';
            fem                     = ElasticProblem(s);
            obj.physicalProblem67 = fem;
        end

        function createElasticProblem78(obj)
            x                       = obj.designVariable;
            s.mesh                  = obj.mesh;
            s.scale                 = 'MACRO';
            s.material              = obj.createMaterial(x);
            s.dim                   = '2D';
            s.boundaryConditions    = obj.createBoundaryConditions(7/9,8/9);
            s.interpolationType     = 'LINEAR';
            s.solverType            = 'REDUCED';
            s.solverMode            = 'DISP';
            s.solverCase = 'rMINRES';
            fem                     = ElasticProblem(s);
            obj.physicalProblem78 = fem;
        end

        function createElasticProblemRight(obj)
            x                        = obj.designVariable;
            s.mesh                   = obj.mesh;
            s.scale                  = 'MACRO';
            s.material               = obj.createMaterial(x);
            s.dim                    = '2D';
            s.boundaryConditions     = obj.createBoundaryConditions(8/9,1);
            s.interpolationType      = 'LINEAR';
            s.solverType             = 'REDUCED';
            s.solverMode             = 'DISP';
            s.solverCase = 'rMINRES';
            fem                      = ElasticProblem(s);
            obj.physicalProblemRight = fem;
        end

        function computeTargetCompliance(obj)
            x     = obj.designVariableTar;
            mat   = obj.createMaterial(x);
            C     = mat.obtainTensor();
            dC    = mat.obtainTensorDerivative();

            %cL     = obj.createComplianceFromConstiutive(obj.physicalProblemLeft);
            cC     = obj.createComplianceFromConstiutive(obj.physicalProblemCenter);
            %cR     = obj.createComplianceFromConstiutive(obj.physicalProblemRight);
            %[jL,~] = cL.computeFunctionAndGradient(C,dC);
            [jC,~] = cC.computeFunctionAndGradient(C,dC);
            %[jR,~] = cR.computeFunctionAndGradient(C,dC);
            obj.targetCompliance = 0.7*jC;
        end

        function c = createComplianceFromConstiutive(obj,physicalProblem)
            s.mesh         = obj.mesh;
            s.stateProblem = physicalProblem;
            c              = ComplianceFromConstiutiveTensor(s);
        end

        function createComplianceLeft(obj)
            x                            = obj.designVariable;
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filter;
            s.complainceFromConstitutive = obj.createComplianceFromConstiutive(obj.physicalProblemLeft);
            s.material                   = obj.createMaterial(x);
            s.value0                     = obj.targetCompliance;
            s.complianceTarget           = 1;
            c                            = ComplianceConstraint(s);
            obj.complianceLeft           = c;
        end

        function createCompliance12(obj)
            x                            = obj.designVariable;
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filter;
            s.complainceFromConstitutive = obj.createComplianceFromConstiutive(obj.physicalProblem12);
            s.material                   = obj.createMaterial(x);
            s.value0                     = obj.targetCompliance;
            s.complianceTarget           = 1;
            c                            = ComplianceConstraint(s);
            obj.compliance12           = c;
        end

        function createCompliance23(obj)
            x                            = obj.designVariable;
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filter;
            s.complainceFromConstitutive = obj.createComplianceFromConstiutive(obj.physicalProblem23);
            s.material                   = obj.createMaterial(x);
            s.value0                     = obj.targetCompliance;
            s.complianceTarget           = 1;
            c                            = ComplianceConstraint(s);
            obj.compliance23           = c;
        end

        function createCompliance34(obj)
            x                            = obj.designVariable;
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filter;
            s.complainceFromConstitutive = obj.createComplianceFromConstiutive(obj.physicalProblem34);
            s.material                   = obj.createMaterial(x);
            s.value0                     = obj.targetCompliance;
            s.complianceTarget           = 1;
            c                            = ComplianceConstraint(s);
            obj.compliance34           = c;
        end

        function createComplianceCenter(obj)
            x                            = obj.designVariable;
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filter;
            s.complainceFromConstitutive = obj.createComplianceFromConstiutive(obj.physicalProblemCenter);
            s.material                   = obj.createMaterial(x);
            s.value0                     = obj.targetCompliance;
            s.complianceTarget           = 1;
            c                            = ComplianceConstraint(s);
            obj.complianceCenter         = c;
        end

        function createCompliance56(obj)
            x                            = obj.designVariable;
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filter;
            s.complainceFromConstitutive = obj.createComplianceFromConstiutive(obj.physicalProblem56);
            s.material                   = obj.createMaterial(x);
            s.value0                     = obj.targetCompliance;
            s.complianceTarget           = 1;
            c                            = ComplianceConstraint(s);
            obj.compliance56           = c;
        end

        function createCompliance67(obj)
            x                            = obj.designVariable;
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filter;
            s.complainceFromConstitutive = obj.createComplianceFromConstiutive(obj.physicalProblem67);
            s.material                   = obj.createMaterial(x);
            s.value0                     = obj.targetCompliance;
            s.complianceTarget           = 1;
            c                            = ComplianceConstraint(s);
            obj.compliance67           = c;
        end

        function createCompliance78(obj)
            x                            = obj.designVariable;
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filter;
            s.complainceFromConstitutive = obj.createComplianceFromConstiutive(obj.physicalProblem78);
            s.material                   = obj.createMaterial(x);
            s.value0                     = obj.targetCompliance;
            s.complianceTarget           = 1;
            c                            = ComplianceConstraint(s);
            obj.compliance78           = c;
        end

        function createComplianceRight(obj)
            x                            = obj.designVariable;
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filter;
            s.complainceFromConstitutive = obj.createComplianceFromConstiutive(obj.physicalProblemRight);
            s.material                   = obj.createMaterial(x);
            s.value0                     = obj.targetCompliance;
            s.complianceTarget           = 1;
            c                            = ComplianceConstraint(s);
            obj.complianceRight          = c;
        end

        function createVolumeFunctional(obj)
            s.mesh         = obj.mesh;
            s.filter       = obj.filter;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            v              = VolumeFunctional(s);
            obj.volume     = v;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.volume;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            h = obj.mesh.computeMinCellSize();
            n = obj.mesh.nnodes;
            M = h^2*sparse(1:n,1:n,ones(n,1));
        end

        function createConstraint(obj)
            switch obj.nLoads
                case 1
                    s.shapeFunctions{1} = obj.complianceCenter;
                case 3
                    s.shapeFunctions{1} = obj.complianceLeft;
                    s.shapeFunctions{2} = obj.complianceCenter;
                    s.shapeFunctions{3} = obj.complianceRight;
                case 9
                    s.shapeFunctions{1} = obj.complianceLeft;
                    s.shapeFunctions{2} = obj.compliance12;
                    s.shapeFunctions{3} = obj.compliance23;
                    s.shapeFunctions{4} = obj.compliance34;
                    s.shapeFunctions{5} = obj.complianceCenter;
                    s.shapeFunctions{6} = obj.compliance56;
                    s.shapeFunctions{7} = obj.compliance67;
                    s.shapeFunctions{8} = obj.compliance78;
                    s.shapeFunctions{9} = obj.complianceRight;
            end
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createDualVariable(obj)
            s.nConstraints   = obj.nLoads;
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
            s.constraintCase = repmat({'INEQUALITY'},[1,obj.nLoads]);
            s.primal         = 'SLERP';
            s.ub             = inf;
            s.lb             = -inf;
            s.rho            = obj.rho;
            opt = OptimizerAugmentedLagrangian(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function m = createMaterial(obj,x)
            f = x.obtainDomainFunction();
            f = obj.filter.compute(f,1);            
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            m = Material.create(s);
        end

        function bc = createBoundaryConditions(obj,t1,t2)
            xMax     = max(obj.mesh.coord(:,1));
            yMax     = max(obj.mesh.coord(:,2));
            isDirX1   = @(coor) abs(coor(:,1))<=0.1*xMax;
            isDirX2   = @(coor) abs(coor(:,1))>=0.9*xMax;
            isDirY   = @(coor) abs(coor(:,2))==0;
            isDir1    = @(coor) isDirX1(coor) & isDirY(coor);
            isDir2    = @(coor) isDirX2(coor) & isDirY(coor);
            isForceX = @(coor) abs(coor(:,1))>=t1*xMax & abs(coor(:,1))<=t2*xMax;
            isForceY = @(coor) abs(coor(:,2))==yMax;
            isForce  = @(coor) isForceX(coor) & isForceY(coor);

            sDir{1}.domain    = @(coor) isDir1(coor);
            sDir{1}.direction = 2;
            sDir{1}.value     = 0;

            sDir{2}.domain    = @(coor) isDir2(coor);
            sDir{2}.direction = [1,2];
            sDir{2}.value     = 0;

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