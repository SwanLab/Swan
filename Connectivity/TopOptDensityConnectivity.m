classdef TopOptDensityConnectivity < handle

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
    end

    methods (Access = public)

        function obj = TopOptDensityConnectivity()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createFilterCompliance();
            obj.createFilterConnectivity();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createEigenValueConstraint();                             
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
            x1      = linspace(0,2,100);
            x2      = linspace(0,1,50);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh.create(s);
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);
            sD.fun      = aFun.project('P1');
            sD.mesh     = obj.mesh;
            sD.type     = 'Density';
            sD.plotting = true;
            dens        = DesignVariable.create(sD);
            obj.designVariable = dens;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filter = f;
        end

        function createFilterCompliance(obj)
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filterComp = f;
%             s.filterType = 'FilterAndProject';
%             s.mesh       = obj.mesh;
%             s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.filterStep = 'PDE';
%             s.beta       = 2.0;
% %             s.eta        = 0.5;
%             f            = Filter.create(s);
%             obj.filterComp = f;
%             s.filterType = 'FilterAdjointAndProject';    
%             f            = Filter.create(s);
%             obj.filterAdjointComp = f;
        end


        function createFilterConnectivity(obj)
            s.filterType = 'FilterAndProject';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            s.filterStep = 'LUMP';
            s.beta       = 100.0;
            s.eta        = 0.05;
%             s.beta = 16.0;
            f            = Filter.create(s);
            obj.filterConnect = f;
            s.filterType = 'FilterAdjointAndProject';    
            f            = Filter.create(s);
            obj.filterAdjointConnect = f;
        end

        function createMaterialInterpolator(obj)
            E0 = 1e-3;
            nu0 = 1/3;
            ndim = obj.mesh.ndim;
            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            E1 = 1;
            nu1 = 1/3;
            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
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
            s.dim                  = '2D';
            m = Material.create(s);
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';
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
            s.filter                      = obj.filterComp;
            if ~isempty(obj.filterAdjointComp)
                s.filterAdjoint               = obj.filterAdjointComp;
            end
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createEigenValueConstraint(obj)                           
            s.mesh              = obj.mesh;
            s.designVariable    = obj.designVariable;
            s.filter            = obj.filterConnect;
            s.filterAdjoint     = obj.filterAdjointConnect;
            s.targetEigenValue = 0.05;      
            s.shift             = 1.0;
            obj.minimumEigenValue = StiffnesEigenModesConstraint(s);
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.shapeFunctions{2} = obj.minimumEigenValue; 
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
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
        end

        function createDualVariable(obj)
            s.nConstraints   = 2;                                        
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
            s.constraintCase{1} = 'EQUALITY';
            s.constraintCase{2} = 'INEQUALITY';                             
            s.ub             = 1;
            s.lb             = 0;
            opt              = OptimizerMMA(s);
            opt.solveProblem();
            obj.optimizer = opt;

%             iterTotal = 0.0;
%             target = 0.06;
%             while target <= 0.15
%                 iterTotal
%                 target
%                 obj.minimumEigenValue.updateTargetEigenvalue(target)
%                 opt              = OptimizerMMA(s);
%                 opt.solveProblem();
%                 target = target + 0.01;
%                 iterTotal = iterTotal + 50;
%             end


%             iterTotal = 400.0;
%             s.maxIter        = 20;
%             beta = 2;
%             while beta <= 64
%                 iterTotal
%                 beta = beta*2;
%                 obj.compliance.updateFilterParams(beta)
%                 opt              = OptimizerMMA(s);
%                 opt.solveProblem();
%                 iterTotal = iterTotal + 20;
%             end
%             
            s.maxIter        = 1000 - iterTotal;
            opt              = OptimizerMMA(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj)
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
