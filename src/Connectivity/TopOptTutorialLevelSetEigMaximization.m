classdef TopOptTutorialLevelSetEigMaximization < handle

    properties (Access = private)
        mesh
        filter
        designVariable
        materialInterpolator
        volume
        cost
        constraint
        dualVariable
        optimizer
        eigenvalue
        dofsNonDesign
        filterAdjoint
        beta
        eta
    end 

    methods (Access = public)
        function obj = TopOptTutorialLevelSetEigMaximization()
            for beta = [1.0, 2.0, 5.0]
                for eta = [0.0, 0.5, 1.0]
                    obj.beta = beta;
                    obj.eta = eta;
                    obj.init()
                    obj.createMesh();
                    obj.createDesignVariable();
                    obj.createFilter();
                    obj.createEigenValue();      
                    obj.createNonDesignableDomain();
                    obj.createVolumeConstraint();
                    obj.createCost();
                    obj.createConstraint();
                    obj.createDualVariable();
                    obj.createOptimizer();
                end
            end
        end
    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            x1      = linspace(0,1,50);
            x2      = linspace(0,1,50);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh.create(s);
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
%             s.filterType = 'LUMP';
%             s.mesh       = obj.mesh;
%             s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
%             f            = Filter.create(s);
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
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.eigenvalue;
            s.weights           = [1.0];
            s.Msmooth           = obj.createMassMatrix();
            s.dofsNonDesign     = obj.dofsNonDesign;
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
            s.dofsNonDesign     = obj.dofsNonDesign;
            obj.constraint      = Constraint(s);
        end

        function createDualVariable(obj)
            s.nConstraints   = 1;
            l                = DualVariable(s);
            obj.dualVariable = l;
        end

        function createOptimizer(obj)
            s.shallPrint     = true;
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.dofsNonDesign  = obj.dofsNonDesign;
            s.GIFname        = 'gif.GIF';
            s.maxIter        = 1000;
            s.tolerance      = 1e-8;
            s.constraintCase{1} = 'EQUALITY';
            s.primal         = 'SLERP';
            s.ub             = inf;
            s.lb             = -inf;
            s.etaNorm        = 0.5; 
            s.etaNormMin     = 0.02;
            s.gJFlowRatio    = 0.2;
            s.etaMax         = 10.0;  
            s.etaMaxMin      = 0.01; 
            s.filter         = obj.filter;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;

%             saveas(figure(1),'1e-35_p2_eta'+string(obj.eta)+'beta'+string(obj.beta)+'.png','png')
            saveas(figure(2),'1e-35_p2_eta'+string(obj.eta)+'beta'+string(obj.beta)+'graficos.png','png')
            writematrix(obj.designVariable.fun.fValues, '1e-35_p2_eta'+string(obj.eta)+'beta'+string(obj.beta)+'.txt')
            f = obj.designVariable.fun.fValues;
            uMesh = obj.designVariable.getUnfittedMesh();
            uMesh.compute(f);
            Fig = figure;
            uMesh.plotStructureInColor('black');
            saveas(figure(Fig), '1e-35_p2_eta'+string(obj.eta)+'beta'+string(obj.beta)+'.png','png')
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
%             isNonDesign =  @(coor) ((abs(coor(:,1)) >= 0.48 & abs(coor(:,1)) <= 0.52) & (abs(coor(:,2)) >= 0.48 & abs(coor(:,2)) <= 0.52)) | ((abs(coor(:,1)) >= 0.25 & abs(coor(:,1)) <= 0.29) & (abs(coor(:,2)) >= 0.25 & abs(coor(:,2)) <= 0.29))  | ((abs(coor(:,1)) >= 0.71 & abs(coor(:,1)) <= 0.75) & (abs(coor(:,2)) >= 0.71 & abs(coor(:,2)) <= 0.75));
            isNonDesign =  @(coor) ((abs(coor(:,1)) >= 0.48 & abs(coor(:,1)) <= 0.52) & (abs(coor(:,2)) >= 0.48 & abs(coor(:,2)) <= 0.52));
            obj.dofsNonDesign = isNonDesign(obj.mesh.coord);
        end


    end
end