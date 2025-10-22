classdef TopOptTutorialLevelSetEigMaximization < handle

    properties (Access = private)
        mesh
        filter
        filterAdjoint
        designVariable
        conductivityInterpolator
        massInterpolator
        boundaryConditions
        eigenvalue
        thermalProblem
        maxTemperature
        volume
        cost
        constraint
        primalUpdater
        optimizer
        dofsNonDesign
        etaNorm
        etaMax
        gJFlowRatio
        perimeter
    end 

    methods (Access = public)
        function obj = TopOptTutorialLevelSetEigMaximization()
            for etaNorm = [0.01] %[0.005] 
                for etaMax = [100.0] %0.15]0.05,0.15,0.5,1,3,4,"cantilever",
                    for gJ = [2.0]
                        obj.etaNorm = etaNorm;
                        obj.etaMax = etaMax;
                        obj.gJFlowRatio = gJ;
                        obj.init()
                        obj.createMesh();
                        obj.createDesignVariable();
                        obj.createFilter();
                        obj.createConductivityInterpolator();
                        obj.createMassInterpolator();
                        obj.createBoundaryConditions();
                        obj.createEigenValue();  
                        obj.createThermalProblem();
                        obj.createMaximumTemperature();
                        obj.createPerimeter();
                        obj.createVolumeConstraint();
                        obj.createCost();
                        obj.createConstraint();
                        obj.createPrimalUpdater();
                        obj.createOptimizer();
                    end
                end
            end
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
%             s.type = 'FourCrossBars';
%             s.xCenter = 0.5;     
%             s.yCenter = 0.5;
%             s.barHalfThickness = 0.05;
%             s.barOffset = 0.30;
%             s.vertHalfLength = 1.0;
%             s.horzHalfLength = 1.0;
            s.type = 'FiveCirclesInclusion';
            s.height = 1.0;
            s.width = 1.0;
            s.radius = sqrt(0.6/(2*pi));
%              s.type = 'FiveCirclesInclusion';
%              s.radius =  sqrt(0.6/(2*pi));
%              s.width = 1.0;
%              s.height = 1.0
%              s.xCoorCenter = 0.0;
%              s.yCoorCenter = 0.0;
%              s.xCoorCenter2 = 1.0;
%              s.yCoorCenter2 = 1.0;
%             s.type = 'Full';
            g      = GeometricalFunction(s);
            lsFun  = g.computeLevelSetFunction(obj.mesh);
            s.fun  = lsFun;
            s.mesh = obj.mesh;
            s.type = 'LevelSet';
            s.plotting = true;
%             s.isFixed.nodes = obj.createNonDesignableDomain();
            ls     = DesignVariable.create(s);
            obj.designVariable = ls;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
%             f.updateEpsilon(2*obj.mesh.computeMinCellSize())
            obj.filter = f;
        end

        function createEigenValue(obj)                           
            s.mesh              = obj.mesh;
            s.designVariable    = obj.designVariable;
            s.filter            = obj.filter;     
            s.filterAdjoint     = obj.filterAdjoint;
            s.conductivityInterpolator = obj.conductivityInterpolator; 
            s.massInterpolator         = obj.massInterpolator; 
            s.isCompl           = true;
            s.dim               = '2D';
            s.boundaryConditions=  obj.boundaryConditions;
            obj.eigenvalue = MaximumEigenValueFunctional(s);
        end

        function createConductivityInterpolator(obj) 
            s.interpolation  = 'SimpAllThermal';
            s.f0   = 1e-2;                                             
            s.f1   = 1;  
            s.dim  = '2D';
            a = MaterialInterpolator.create(s);
            obj.conductivityInterpolator = a;            
        end 

        function createMassInterpolator(obj)
            s.interpolation  = 'SIMPThermal';                              
            s.f0   = 1e-2;
            s.f1   = 1;
            s.pExp = 1;
            a = MaterialInterpolator.create(s);
            obj.massInterpolator = a;            
        end      

       function createThermalProblem(obj)
            s.mesh = obj.mesh;
            s.conductivity = obj.conductivityInterpolator; 
            s.mass = obj.massInterpolator; 
            Q = LagrangianFunction.create(obj.mesh,1,'P1');
            fValues = ones(Q.nDofs,1);
            Q.setFValues(fValues);
            s.source = Q;  
            s.dim = '2D';
            s.boundaryConditions =  obj.boundaryConditions;
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';
            fem = ThermalProblem(s); 
            obj.thermalProblem = fem;
       end

       function createMaximumTemperature(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filter;
            s.stateProblem                = obj.thermalProblem;
            s.conductivity                =  obj.conductivityInterpolator; 
            s.mass                        =  obj.massInterpolator; 
            c = MaximumTemperatureFunctional(s);  
            obj.maxTemperature = c;
        end

        function uMesh = createBaseDomain(obj)
            sG.type          = 'Full';
            g                = GeometricalFunction(sG);
            lsFun            = g.computeLevelSetFunction(obj.mesh);
            levelSet         = lsFun.fValues;
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh            = UnfittedMesh(s);
            uMesh.compute(levelSet);
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.test = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;
            s.uMesh = obj.createBaseDomain();
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createPerimeter(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filter;
            p = SimplePerimeterFunctional(s);
            obj.perimeter = p;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.eigenvalue;
            s.shapeFunctions{2} = obj.perimeter;
            s.weights           = [1.0,10.0];
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

        function createPrimalUpdater(obj)
            s.mesh = obj.mesh;
            obj.primalUpdater = SLERP(s);
        end

        function createOptimizer(obj)
            s.GIFname        = 'lsEigMax';
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 1500;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.primalUpdater  = obj.primalUpdater;
            s.etaNorm        = obj.etaNorm;
            s.etaNormMin     = obj.etaNorm;
            s.gJFlowRatio    = obj.gJFlowRatio;
            s.etaMax         = obj.etaMax;
            s.etaMaxMin      = 0.02;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;

            saveas(figure(1),'design'+string(obj.etaNorm)+'-'+string(obj.etaMax)+'-'+string(obj.gJFlowRatio)+'-'+'.png','png')
            saveas(figure(2),'graficos'+string(obj.etaNorm)+'-'+string(obj.etaMax)+'-'+string(obj.gJFlowRatio)+'-'+'.png','png')
            writematrix(obj.designVariable.fun.fValues, 'fvalues.txt')
            f = obj.designVariable.fun.fValues;
            uMesh = obj.designVariable.getUnfittedMesh();
            uMesh.compute(f);
            Fig = figure;
            uMesh.plotStructureInColor('black');
            saveas(figure(Fig), 'designBlack.png','png')
        end

        function createBoundaryConditions(obj)
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isLeft  = @(coor) abs(coor(:,1))==xMin;
            isRight = @(coor) abs(coor(:,1))==xMax;
            isDown = @(coor) abs(coor(:,2))==yMin;
            isUp = @(coor) abs(coor(:,2))== yMax;
            isDir   = @(coor) isLeft(coor) | isRight(coor) | isUp(coor) | isDown(coor);  
%             isDir = @(coor) isRight(coor) | isUp(coor);
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
            obj.boundaryConditions = BoundaryConditions(s);  
        end

        function [dofsNonDesign] = createNonDesignableDomain(obj)
%             isNonDesign =  @(coor) ((abs(coor(:,1)) >= 0.48 & abs(coor(:,1)) <= 0.52) & (abs(coor(:,2)) >= 0.48 & abs(coor(:,2)) <= 0.52)) | ((abs(coor(:,1)) >= 0.25 & abs(coor(:,1)) <= 0.29) & (abs(coor(:,2)) >= 0.25 & abs(coor(:,2)) <= 0.29))  | ((abs(coor(:,1)) >= 0.71 & abs(coor(:,1)) <= 0.75) & (abs(coor(:,2)) >= 0.71 & abs(coor(:,2)) <= 0.75));
%             isNonDesign =  @(coor) ((abs(coor(:,1)) >= 0.48 & abs(coor(:,1)) <= 0.52) & (abs(coor(:,2)) >= 0.48 & abs(coor(:,2)) <= 0.52));
            r = 0.08;
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
%             centers = [0.1,0.1;
%                        0.1,0.9;
%                        0.9,0.1;
%                        0.9,0.9];
            centers = [0.5 0.5];
            
            isNonDesign = @(coor) any( (coor(:,1) - centers(:,1)').^2 + ...
                                       (coor(:,2) - centers(:,2)').^2 <= r^2, 2);

            dofsNonDesign = isNonDesign(obj.mesh.coord);
        end


    end
end