classdef MeshSymmetrizerTest < testShowingError

    properties (Access = protected)
        tol = 1e-14;
    end

    properties (Access = private)
        mesh
        symmetricMesh
    end

    properties (Access = private)
        alpha
        nx
        ny
        x1Max
        x1Min
        x2Max
        x2Min
    end

    methods (Access = public)

        function obj = MeshSymmetrizerTest()
            obj.init();
            obj.createMesh();
            obj.createSymmetricMesh();
            obj.plotMeshes();
        end

    end

    methods (Access = private)

        function init(obj)
            obj.nx = 20;
            obj.ny = 30;
            obj.x1Max = 3;
            obj.x1Min = 0;
            obj.x2Max = 2;
            obj.x2Min = 0;
            obj.alpha = 0.3;
        end

        function createMesh(obj)
            [x1,x2] = obj.createCoord();
            x3 = zeros(size(x2));
            [F,V] = mesh2tri(x1,x2,x3,'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            m = Mesh.create(s);
            obj.mesh = m;
        end

        function [x1,x2] = createCoord(obj)
            [x1,x2] = obj.createHorizontalVerticalCoords();
            [x1,x2] = obj.rotateCoord(x1,x2);
        end

        function [x1v,x2v] = createHorizontalVerticalCoords(obj)
            x1 = linspace(obj.x1Min,obj.x1Max,obj.nx);
            x2 = linspace(obj.x2Min,obj.x2Max,obj.ny);
            [x1v,x2v] = meshgrid(x1,x2);
        end

        function [x1r,x2r] = rotateCoord(obj,x1,x2)
            ca = cos(obj.alpha);
            sa = sin(obj.alpha);
            x1r = ca*x1 - sa*x2;
            x2r = sa*x1 + ca*x2;
        end

        function createSymmetricMesh(obj)
            s.mesh = obj.mesh;
            s.symmetricLine.vector = [cos(obj.alpha);sin(obj.alpha)];
            s.symmetricLine.point = [0;0];
            mS = Symmetrizer(s);
            m = mS.computeSymmetricMesh();
            obj.symmetricMesh = m;
        end

        function plotMeshes(obj)
            figure()
            obj.mesh.plot();
            figure()
            obj.symmetricMesh.plot();
        end

    end

    methods (Access = protected)

        function computeError(obj)
            obj.error = 0;
        end

    end

end