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
        primalUpdater
    end 

    methods (Access = public)
        function obj = TopOptTutorialLevelSetEigMaximization()
                    obj.init()
                    obj.createMesh();
                    obj.createDesignVariable();
%                     obj.createFilter();
                    obj.createEigenValue();      
                    obj.createNonDesignableDomain();
                    obj.createVolumeConstraint();
                    obj.createCost();
                    obj.createConstraint();
                    obj.createPrimalUpdater();
                    obj.createOptimizer();
%                 end
%             end
        end
    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            x1      = linspace(0,1.0,50);
            x2      = linspace(0,1.0,50);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh.create(s);
        end

        function createDesignVariable(obj)
%             s.type        = 'ThreeRectangles';
%             s.xSide1       = 0.3;
%             s.ySide1       = 0.3;
%             s.xCoorCenter1 = 0.35;
%             s.yCoorCenter1 = 0.5;
%             s.xSide2       = 0.1;
%             s.ySide2       = 0.1;
%             s.xCoorCenter2 = 0.6;
%             s.yCoorCenter2 = 0.5;
%             s.xSide3       = 1.0;
%             s.ySide3       = 0.1;
%             s.xCoorCenter3 = 0.5;
%             s.yCoorCenter3 = 0.2; 
%             s.type = 'TwoCircles';
%             s.radius = 0.05;
%             s.xCoorCenter = 0.1;
%             s.yCoorCenter = 0.1;
%             s.xCoorCenter2 = 0.3;
%             s.yCoorCenter2 = 0.3;
            s.type = 'Circle';
            s.radius = sqrt(0.4/pi);
            s.xCoorCenter = 0.5;
            s.yCoorCenter = 0.5;      
%             s.type = 'Full';
            g      = GeometricalFunction(s);
            lsFun  = g.computeLevelSetFunction(obj.mesh);
            s.fun  = lsFun;
            s.mesh = obj.mesh;
            s.type = 'LevelSet';
            s.plotting = true;
            s.isFixed.nodes = obj.createNonDesignableDomain();
            ls     = DesignVariable.create(s);
%             ls.fun.setFValues(importdata('fvalues2.txt'))
            obj.designVariable = ls;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filter = f;
        end

        function createEigenValue(obj)                           
            s.mesh              = obj.mesh;
            s.designVariable    = obj.designVariable;
            s.filter            = obj.filter;     
            s.filterAdjoint     = obj.filterAdjoint;
            s.boundaryConditions= obj.createEigenvalueBoundaryConditions();
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
            n = obj.mesh.nnodes;
            h = obj.mesh.computeMinCellSize();
            M = h^2*sparse(1:n,1:n,ones(1,n),n,n);
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();
            s.dofsNonDesign     = obj.dofsNonDesign;
            obj.constraint      = Constraint(s);
        end

        function createPrimalUpdater(obj)
            s.mesh = obj.mesh;
            obj.primalUpdater = SLERP(s);
        end

        function createOptimizer(obj)
            s.shallPrint     = true;
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.primalUpdater  = obj.primalUpdater;
            s.GIFname        = 'gif.GIF';
            s.maxIter        = 1000;
            s.tolerance      = 1e-8;
            s.constraintCase{1} = 'EQUALITY';
            s.primal         = 'SLERP';
            s.ub             = inf;
            s.lb             = -inf;
            s.etaNorm        = 0.02; %0.02
            s.etaNormMin     = 0.02;
            s.gJFlowRatio    = 1.0;  %0.2   2.0; 60.0
            s.etaMax         = 1.0; %1.0; %60.0;    % 1 - 5.0 5.0
            s.etaMaxMin      = 0.01; %0.01; %10.0; 
            s.filter         = obj.filter;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;

            saveas(figure(1),'design.png','png')
            saveas(figure(2),'graficos.png','png')
            writematrix(obj.designVariable.fun.fValues, 'fvalues.txt')
            f = obj.designVariable.fun.fValues;
            uMesh = obj.designVariable.getUnfittedMesh();
            uMesh.compute(f);
            Fig = figure;
            uMesh.plotStructureInColor('black');
            saveas(figure(Fig), 'designBlack.png','png')
        end

        function  bc = createEigenvalueBoundaryConditions(obj)
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isLeft  = @(coor) abs(coor(:,1))==xMin;
            isRight = @(coor) abs(coor(:,1))==xMax;
            isDown = @(coor) abs(coor(:,2))==yMin;
            isUp = @(coor) abs(coor(:,2))== yMax;
            isDir   = @(coor) isLeft(coor) | isRight(coor) | isUp(coor) | isDown(coor);  
%             isDir   = @(coor)   isRight(coor) | isUp(coor);  
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

        function [dofsNonDesign] = createNonDesignableDomain(obj)
%             isNonDesign =  @(coor) ((abs(coor(:,1)) >= 0.48 & abs(coor(:,1)) <= 0.52) & (abs(coor(:,2)) >= 0.48 & abs(coor(:,2)) <= 0.52)) | ((abs(coor(:,1)) >= 0.25 & abs(coor(:,1)) <= 0.29) & (abs(coor(:,2)) >= 0.25 & abs(coor(:,2)) <= 0.29))  | ((abs(coor(:,1)) >= 0.71 & abs(coor(:,1)) <= 0.75) & (abs(coor(:,2)) >= 0.71 & abs(coor(:,2)) <= 0.75));
%             isNonDesign =  @(coor) ((abs(coor(:,1)) >= 0.48 & abs(coor(:,1)) <= 0.52) & (abs(coor(:,2)) >= 0.48 & abs(coor(:,2)) <= 0.52));
            r = 0.1;
%             x0 = 0.5; y0 = 0.5;
%             isNonDesign =  @(coor) (((coor(:,1)-x0).^2+(coor(:,2)-y0).^2-r^2 <= 0));
%             r = 0.05;
%             x0 = 0.1; y0 = 0.1;
%             x02 = 0.3; y02 = 0.3;
%             isNonDesign =  @(coor) (((coor(:,1)-x0).^2+(coor(:,2)-y0).^2-r^2 <= 0) | ((coor(:,1)-x02).^2+(coor(:,2)-y02).^2-r^2 <= 0));

%             centers = [0.5 0.5;
%                        0.25 0.25;
%                        0.25 0.75;
%                        0.75 0.25;
%                        0.75 0.75];

            centers = [0.5 0.5];
            
            isNonDesign = @(coor) any( (coor(:,1) - centers(:,1)').^2 + ...
                                       (coor(:,2) - centers(:,2)').^2 <= r^2, 2);

            dofsNonDesign = isNonDesign(obj.mesh.coord);
        end


    end
end