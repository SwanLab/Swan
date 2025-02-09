classdef MeshComponentCounterTest < testShowingError

    properties (Access = protected)
        tol = 1e-14;
    end

    properties (Access = private)

    end

    properties (Access = private)
        backgroundMesh
        uMesh
    end

    methods (Access = public)

        function obj = MeshComponentCounterTest()
            obj.init();
            obj.createBackgroundMesh();
            obj.createUnfittedMesh();
            obj.countMeshesComponents();
        end

    end

    methods (Access = private)

        function init(obj)

        end

        function createBackgroundMesh(obj)
            x1 = linspace(0,1,3);
            x2 = linspace(0,1,3);
            [x1,x2] = meshgrid(x1,x2);
            x3 = zeros(size(x2));
            [F,V] = mesh2tri(x1,x2,x3,'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            m = Mesh.create(s);
            obj.backgroundMesh = m;
        end

        function createUnfittedMesh(obj)
            ls = obj.createLevelSet();
            s.boundaryMesh   = obj.createBoundaryMesh();
            s.backgroundMesh = obj.backgroundMesh;
            obj.uMesh = UnfittedMesh(s);
            obj.uMesh.compute(ls);
            %  obj.uMesh.plotStructureInColor('black')
        end

        function ls = createLevelSet(obj)
            s.coord = obj.backgroundMesh.coord;
            s.type   = 'rectangleInclusion';
            s.widthV = 1.1;
            s.widthH = 0.3;
            s.ndim   = 2;
            levelSet = LevelSetCreator.create(s);
            ls = levelSet.getValue();
        end

        function m = createBoundaryMesh(obj)
            sB.backgroundMesh = obj.backgroundMesh;
            sB.dimension = 1:3;
            sB.type = 'FromReactangularBox';
            bMc = BoundaryMeshCreator.create(sB);
            m  = bMc.create();
        end

        function countMeshesComponents(obj)
            s.unfittedMesh = obj.uMesh;
            sp = UnfittedMeshSplitter(s);
            sp.split();
            sp.plot();
        end

    end

    methods (Access = protected)

        function computeError(obj)
            obj.error = 0;
        end

    end



end