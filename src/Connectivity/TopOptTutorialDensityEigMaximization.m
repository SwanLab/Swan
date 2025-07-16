classdef TopOptTutorialDensityEigMaximization < handle

    properties (Access = private)
        mesh
        filter
        filterAdjoint
        designVariable
        materialInterpolator
        volume
        cost
        constraint
        dualVariable
        optimizer
        eigenvalue
        dofsNonDesign
        beta
        eta
    end 

    methods (Access = public)
        function obj = TopOptTutorialDensityEigMaximization()
            obj.beta = 5.0;
            obj.eta = 0.5;
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createNonDesignableDomain();
            obj.createEigenValue();                             
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
            x1      = linspace(0,1,100);
            x2      = linspace(0,1,100);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh.create(s);
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) 1.0*ones(size(x(1,:,:)));
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
%             s.filterType = 'LUMP';
%             s.mesh  = obj.mesh;
%             s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
%             f = Filter.create(s);
%             obj.filter = f;

            s.filterType = 'FilterAndProject';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            s.filterStep = 'LUMP';
            s.beta       = obj.beta;
            s.eta        = obj.eta;
            obj.filter = Filter.create(s);

            s.filterType = 'FilterAdjointAndProject';   
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            s.filterStep = 'LUMP';
            s.beta       = obj.beta;
            s.eta        = obj.eta;
            obj.filterAdjoint = Filter.create(s);

        end

        function createEigenValue(obj)                           
            s.mesh              = obj.mesh;
            s.designVariable    = obj.designVariable;
            s.filter            = obj.filter;
            s.filterAdjoint     = obj.filterAdjoint;
            s.boundaryConditions= obj.createEigenvalueBoundaryConditions();
            s.shift             = 0.0;
            obj.eigenvalue = MaximumEigenValueFunctional(s);
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;
%             s.dofsNonDesign = obj.dofsNonDesign;
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.eigenvalue;
            s.weights           = [1.0];
%             s.dofsNonDesign     = obj.dofsNonDesign;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            s.test  = LagrangianFunction.create(obj.mesh,1,'P1');
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh  = obj.mesh;
            s.type  = 'MassMatrix';
            LHS = LHSIntegrator.create(s);
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
%             s.dofsNonDesign  = obj.dofsNonDesign;
            s.maxIter        = 3000;
            s.tolerance      = 1e-8;
            s.constraintCase{1} = 'EQUALITY';
            s.ub             = 1;
            s.lb             = 0;
            opt              = OptimizerMMA(s);
            opt.solveProblem();
            obj.optimizer = opt;
%             saveas(figure(1),'eigMaxDesign.png','png')
%             saveas(figure(2),'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+'graficos.png','png')
%             writematrix(obj.designVariable.fun.fValues,'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+'.txt')
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

        function createNonDesignableDomain(obj)
            isNonDesign =  @(coor) ((abs(coor(:,1)) >= 0.48 & abs(coor(:,1)) <= 0.52) & (abs(coor(:,2)) >= 0.48 & abs(coor(:,2)) <= 0.52)) | ((abs(coor(:,1)) >= 0.25 & abs(coor(:,1)) <= 0.29) & (abs(coor(:,2)) >= 0.25 & abs(coor(:,2)) <= 0.29))  | ((abs(coor(:,1)) >= 0.75 & abs(coor(:,1)) <= 0.79) & (abs(coor(:,2)) >= 0.75 & abs(coor(:,2)) <= 0.79));
%             isNonDesign =  @(coor) ((abs(coor(:,1)) >= 0.48 & abs(coor(:,1)) <= 0.52) & (abs(coor(:,2)) >= 0.48 & abs(coor(:,2)) <= 0.52));
            obj.dofsNonDesign = isNonDesign(obj.mesh.coord);
        end


    end
end