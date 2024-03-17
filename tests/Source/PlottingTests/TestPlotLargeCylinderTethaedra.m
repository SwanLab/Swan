classdef TestPlotLargeCylinderTethaedra < handle

    properties (Access = private)
        backgroundMesh
        boundaryMesh
        levelSet
        unfittedMesh
    end

    methods (Access = public)

        function obj = TestPlotLargeCylinderTethaedra()
            obj.createBackgroundMesh();
            obj.createBoundaryMesh();
            obj.createLevelSet();
            obj.createUnfittedMesh();
            obj.plotUnfittedMesh();
        end

    end

    methods (Access = private)

        function createBackgroundMesh(obj)
            x = linspace(0,1,10);
            y = linspace(0,1,10);
            z = linspace(0,2,20);
            [X,Y,Z] = meshgrid(x,y,z);
            coord  = [X(:) Y(:) Z(:)];
            d = delaunayTriangulation(coord);
            s.connec = d.ConnectivityList;
            s.coord  = coord;
            obj.backgroundMesh = Mesh.create(s);
        end

        function createBoundaryMesh(obj)
            s.backgroundMesh = obj.backgroundMesh;
            s.dimension = 1:obj.backgroundMesh.ndim;
            bC = BoundaryMeshCreatorFromRectangularBox(s);
            bM = bC.create();
            obj.boundaryMesh = bM;
        end

        function createLevelSet(obj)
            mesh = obj.backgroundMesh;
            s.type = 'Cylinder';
            s.radius = 1.1;
            s.xCoorCenter = 0.5;
            s.yCoorCenter = 0.5;
            fun = GeometricalFunction(s);
            lsfun = fun.computeLevelSetFunction(mesh);
            obj.levelSet = lsfun.fValues;
        end

        function createUnfittedMesh(obj)
            s.backgroundMesh = obj.backgroundMesh;
            s.boundaryMesh   = obj.boundaryMesh;
            uM = UnfittedMesh(s);
            uM.compute(obj.levelSet);
            obj.unfittedMesh = uM;
        end

        function plotUnfittedMesh(obj)
            figure();
            obj.unfittedMesh.plot();
            view([1 1 1])
        end

    end
end

