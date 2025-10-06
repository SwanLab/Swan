classdef ConnectivityComputer < handle

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
        monitoringEigenModes
    end

    properties (Access = private)
        
    end

    methods (Access = public)

        function obj = ConnectivityComputer()
%             lambda1 = []
%             h = 1/100;
% %             for radius = 0:h:0.6 %0.5
%             radius_all = 0.0:h:0.6
%             for radius = radius_all %0.5
%                 radius
            obj.init();
            obj.createMesh();
            obj.createLevelSet(); %radius
            obj.createFilter();
            obj.createCharacteristicFunction();
            obj.createDesignVariable();  
            lambda = obj.computeEigenValueFunctional();
%             lambda1 = [lambda1, lambda]
%             end
%             figure
%             semilogy(radius_all,lambda1,'-')
%             xlabel('Radius')
%             ylabel('First Eigenvalue')
%             grid on
% %             save('PDEFP.mat','lambda1')
%             save('LUMP.mat','lambda1')
% %             save('PDE.mat','lambda1')
        end

    end

    methods (Access = private)

        function init(obj)
            close all
        end

        function createMesh(obj)
            x1 = linspace(0,1,100);
            x2 = linspace(0,1,100);
            [xv,yv] = meshgrid(x1,x2);
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            m = Mesh.create(s);
            obj.mesh = m;
        end

        function createLevelSet(obj) %,radius
%             s.type        = 'CircleInclusion';
%             s.xCoorCenter       = 0.5;
%             s.yCoorCenter       = 0.5;
%             s.radius      = radius;% 
% % 
%             s.type        = 'RingWithHorizontalCrack'; %'RingSDF';
%             s.innerRadius = 0.1;
%             s.outerRadius = 0.2;
%             s.xCoorCenter = 0.5;
%             s.yCoorCenter = 0.5;
%             s.crackWidth = 0.01;
% 
%             s.type        = 'RingSDF';
%             s.innerRadius = 0.1;
%             s.outerRadius = 0.15;
%             s.xCoorCenter = 0.5;
%             s.yCoorCenter = 0.5;

            s.type        = 'ThreeRectanglesInclusion';
            s.xSide1       = 0.3;
            s.ySide1       = 0.3;
            s.xCoorCenter1 = 0.35;
            s.yCoorCenter1 = 0.5;
            s.xSide2       = 0.1;
            s.ySide2       = 0.1;
            s.xCoorCenter2 = 0.6;
            s.yCoorCenter2 = 0.5;
            s.xSide3       = 1.0;
            s.ySide3       = 0.1;
            s.xCoorCenter3 = 0.5;
            s.yCoorCenter3 = 0.2; 

            g             = GeometricalFunction(s);
            phi           = g.computeLevelSetFunction(obj.mesh);
            obj.levelSet = phi;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;

%             s.filterType = 'FilterAndProject';
% %             s.filterType = 'CloseOperator'; %'FilterAndProject';
%             s.mesh       = obj.mesh;
%             s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.filterStep = 'LUMP';
%             s.beta       = 100.0; % 1.0;
%             s.eta        = 0.5;
%             obj.filter = Filter.create(s);
        end        

        function createCharacteristicFunction(obj)
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh;
            uMesh            = UnfittedMesh(s);
            uMesh.compute(obj.levelSet.fValues);
            obj.characteristicFunction  = CharacteristicFunction.create(uMesh);
        end

        function createDesignVariable(obj)
            s.fun  = obj.filter.compute(obj.characteristicFunction,3);
            s.mesh = obj.mesh;
            s.type = 'Density';
            s.plotting = true;
            dens    = DesignVariable.create(s);
            dens.plot();
            obj.designVariable = dens;
        end

        function [lambda] = computeEigenValueFunctional(obj)
            s.mesh = obj.mesh;
            s.designVariable = obj.designVariable;
            s.filter = obj.filter;
            s.boundaryConditions = obj.createEigenvalueBoundaryConditions();
            s.eigenModes = StiffnessEigenModesComputer(s);
            s.isCompl  = true;
            mE = MinimumEigenValueFunctional(s);
            [lambda, dlambda] = mE.computeFunctionAndGradient(obj.designVariable);  
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