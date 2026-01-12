classdef ConnectivityComputerTwoPhaseTest < handle

    properties (Access = public)

    end

    properties (Access = private)
        mesh
        levelSet
        characteristicFunction
        designVariable
        materialInterpolation
        filter
        filterConnect
        filterAdjointConnect
        eigVs
        eigFs
    end

    properties (Access = private)
        
    end

    methods (Access = public)

        function obj = ConnectivityComputerTwoPhaseTest()
            obj.init();
            obj.createMesh();
            obj.createLevelSet();
            obj.createCharacteristicFunction();
            obj.createFilter();
            obj.createDesignVariable(); 
            obj.createFilterConnectivity();
            obj.computeEigenvaluesWithDifferentParameters();
%             obj.createPlot();
        end
    end

    methods (Access = private)

        function init(obj)
            close all
        end

        function createMesh(obj)
            x1 = linspace(0,1.0,50);
            x2 = linspace(0,1.0,50);
            [xv,yv] = meshgrid(x1,x2);
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            m = Mesh.create(s);
            obj.mesh = m;
        end

        function createLevelSet(obj)
%             s.type        = 'RectangleInclusion';
%             s.xSide       = 0.5;
%             s.ySide       = 0.5;
%             s.xCoorCenter = 0.5;
%             s.yCoorCenter = 0.5;
            s.type        = 'Empty';
%             s.type        = 'ThreeRectanglesInclusion';
%             s.xSide1       = 0.3;
%             s.ySide1       = 0.3;
%             s.xCoorCenter1 = 0.4;
%             s.yCoorCenter1 = 0.5;
%             s.xSide2       = 0.05;
%             s.ySide2       = 0.05;
%             s.xCoorCenter2 = 0.6;
%             s.yCoorCenter2 = 0.5;
%             s.xSide3       = 1.0;
%             s.ySide3       = 0.1;
%             s.xCoorCenter3 = 0.5;
%             s.yCoorCenter3 = 0.2; 
            g             = GeometricalFunction(s);
            phi           = g.computeLevelSetFunction(obj.mesh);
            obj.levelSet = phi;
            obj.levelSet.setFValues(importdata('maxCaseC.txt'))
        end

        function createCharacteristicFunction(obj)
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh;
            uMesh            = UnfittedMesh(s);
            uMesh.compute(obj.levelSet.fValues);
            obj.characteristicFunction  = CharacteristicFunction.create(uMesh);
        end

        function createDesignVariable(obj)
%             s.fun  = obj.filter.compute(obj.characteristicFunction,3);
%             s.fValues = round(obj.characteristicFunction.project('P1').fValues);
%             s.mesh    = obj.mesh;
%             s.order   = 'P1';
%             s.fun = LagrangianFunction(s);

%             cant = LagrangianFunction.create(obj.mesh,1,'P1');
%             cant.setFValues(importdata('cantilever.txt'));
%             cant.fValues = importdata('optBridgeConnect.txt');
%             s.fun = cant;    

%             s.mesh = obj.mesh;
%             s.type = 'Density';
%             s.plotting = true;
%             dens    = DesignVariable.create(s);
%             dens.plot();
%             obj.designVariable = dens;

            s.fun  = obj.levelSet;
            s.mesh = obj.mesh;
            s.type = 'LevelSet';
            s.plotting = true;
            ls     = DesignVariable.create(s);
            obj.designVariable = ls;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end       

        function createFilterConnectivity(obj)
%             s.filterType = 'PDE';
%             s.mesh  = obj.mesh;
%             s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
%             f = Filter.create(s);
% %             f.updateEpsilon(3.0*obj   .mesh.computeMeanCellSize());
%             obj.filterConnect = f;
% 
%             s.filterType = 'FilterAndProject';
%             s.mesh       = obj.mesh;
%             s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.filterStep = 'LUMP';
%             s.beta       = 10.0;
%             s.eta        = 0.5;
%             f            = Filter.create(s);
%             obj.filterConnect = f;
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filterConnect = f;
        end        

        function  [lambdas, phis] = computeEigenValueFunctional(obj, n, epsilon, p)
            eigen = obj.computeEigenValueProblem(epsilon, p);
            s.eigenModes = eigen;
            s.designVariable = obj.designVariable;
            s.filter = obj.filterConnect;
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.createEigenvalueBoundaryConditions();
%             mE = MinimumEigenValueFunctional(s);
            mE = MaximumEigenValueFunctional(s);
            [lambdas, phis] = mE.computeEigenModes(obj.designVariable, n);
        end

        function eigen = computeEigenValueProblem(obj,epsilon, p)
            s.mesh  = obj.mesh;
            s.epsilon = epsilon;
            s.p       = p;
            s.boundaryConditions = obj.createEigenvalueBoundaryConditions();
            eigen   = StiffnessEigenModesComputer(s);
        end

        function createPlot(obj)
            analytical = [0.0, 4*pi^2,  4*pi^2,  8*pi^2,  16*pi^2,  16*pi^2 ];
            error = abs(analytical- obj.eigVs);
            error(:,2:6) = error(:,2:6)./analytical(2:6);
            figure;
            loglog([1e-3, 1e-5, 1e-7, 1e-9, 1e-12], error,'o-');
            grid; ylim([1e-12,1e0])
            legend({'\lambda_1', '\lambda_2','\lambda_3','\lambda_4','\lambda_5','\lambda_6'},'Location','southeast');
            xlabel('\alpha_0 = \beta_0')
            ylabel('Error')
%             figure;
%             semilogy([1, 2, 4, 8], error,'o-');
%             legend({'\lambda_1', '\lambda_2','\lambda_3','\lambda_4','\lambda_5','\lambda_6'},'Location','southeast');
%             xlabel('p');
%             ylabel('Error');
%             grid; ylim([1e-6,1e3])
%             title('\alpha_0 = \beta_0 = 1e-5, LUMP')
        end

        function computeEigenvaluesWithDifferentParameters(obj)
            n = 12;
            j = 1;
            eigVs= [];
            eigFs= []; 
            for epsilon = [1e-5]%[1e-3, 1e-5, 1e-7, 1e-9, 1e-12]
%             for p = [1, 2, 4, 8]
                [lambda, phi] = obj.computeEigenValueFunctional(n, epsilon, 8);
                eigVs = [eigVs; lambda'];
                eigFs = [eigFs; phi];
                eigVs'
                j = j + 1;
            end
            obj.eigVs = eigVs;
            obj.eigFs = eigFs;
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

    end

end
