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
    end

    properties (Access = private)
        
    end

    methods (Access = public)

        function obj = ConnectivityComputer()
            obj.init();
            obj.createMesh();
            obj.createLevelSet();
            obj.createFilter();
            obj.createCharacteristicFunction();
            obj.createDesignVariable();  
            obj.computeEigenValueFunctional()
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

        function createLevelSet(obj)
%             s.type        = 'CircleInclusion';
%             s.radius      = 0.4;
%             s.xCoorCenter = 0.5;
%             s.yCoorCenter = 0.5;
%             s.type = 'Empty'    
            s.type        = 'RectangleInclusion';
            s.xSide       = 0.5;
            s.ySide       = 0.5;
            s.xCoorCenter = 0.5;
            s.yCoorCenter = 0.5;
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

        function computeEigenValueFunctional(obj)
            eigen = obj.computeEigenValueProblem();
            s.eigenModes = eigen;
            s.designVariable = obj.designVariable;
            s.mesh = obj.mesh;
            mE = MinimumEigenValueFunctional(s);
            [lambda, dlambda] = mE.computeFunctionAndGradient(obj.designVariable);
            dlambda = obj.filter.compute(dlambda,2);
            dlambda.plot()
        end

        function eigen = computeEigenValueProblem(obj)
            s.mesh  = obj.mesh;
            s.shift = 1.0;
            eigen   = StiffnessEigenModesComputer(s);
        end



    end

end
