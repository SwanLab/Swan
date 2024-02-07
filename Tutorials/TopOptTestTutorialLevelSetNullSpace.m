classdef TopOptTestTutorialLevelSetNullSpace < handle

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
        primalUpdater
        optimizer
    end

    methods (Access = public)

        function obj = TopOptTestTutorialLevelSetNullSpace()
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
            obj.createPrimalUpdater();
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj)

        end

        function createMesh(obj)
            %UnitMesh better
            x1      = linspace(0,1,50);
            x2      = linspace(0,1,50);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh(s);            
        end

        function createDesignVariable(obj)
            s.type = 'Full';
            g      = GeometricalFunction(s);
            lsFun  = g.computeLevelSetFunction(obj.mesh);           
            s.fun  = lsFun;
            s.mesh = obj.mesh;                        
            s.type = 'LevelSet';
            ls     = DesignVariable.create(s);   
            obj.designVariable = ls;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = P1Function.create(obj.mesh,1);
            f = Filter.create(s);
            obj.filter = f;
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

        function createElasticProblem(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = f.project('P1');
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createInterpolatedMaterial(f);
            s.dim = '2D';
            s.bc = obj.createBoundaryConditions();
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
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filter;
            s.complainceFromConstitutive = obj.createComplianceFromConstiutive();
            s.materialInterpolator       = obj.materialInterpolator;
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.volumeTarget = 0.4;
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function V = createVolumeFunctional(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            V        = VolumeFunctional(s);
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.weights           = 1;
            obj.cost            = Cost(s);
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            obj.constraint      = Constraint(s);
        end

        function createDualVariable(obj)
            s.nConstraints   = 1;
            l                = DualVariable(s);
            obj.dualVariable = l;
        end

        function m = createMonitoring(obj,isCreated)
            switch isCreated
                case true
                    s.type = 'OptimizationProblem';
                case false
                    s.type = 'Null';
            end
            s.cost             = obj.cost;
            s.constraint       = obj.constraint;
            s.designVariable   = obj.designVariable;
            s.dualVariable     = obj.dualVariable;
            s.functionals{1}   = obj.createVolumeFunctional();
            s.optimizationType = 'NullSpace';
            s.primalUpdater    = obj.primalUpdater;
            s.isConstrained    = true;
            s.maxNColumns      = 5;
            m                  = Monitoring.create(s);
        end

        function createPrimalUpdater(obj)
            s.primal          = 'SLERP';
            s.designVariable  = obj.designVariable;
            f                 = PrimalUpdaterFactory();
            obj.primalUpdater = f.create(s);
        end

        function createOptimizer(obj)
            s.monitoring     = obj.createMonitoring(true);
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 100;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.volumeTarget   = 0.4;
            s.primal         = 'SLERP';
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function mat = createInterpolatedMaterial(obj,dens)
            mI   = obj.materialInterpolator;
            mat  = mI.computeConsitutiveTensor(dens);
        end
        
        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.3*yMax & abs(coor(:,2))<=0.7*yMax);

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
            bc.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = PointLoad(obj.mesh, sPL{i});
                pointloadFun = [pointloadFun, pl];
            end
            bc.pointloadFun = pointloadFun;

            bc.periodicFun  = [];
        end
    end
end