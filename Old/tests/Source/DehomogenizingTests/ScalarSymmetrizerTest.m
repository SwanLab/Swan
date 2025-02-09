classdef ScalarSymmetrizerTest < testShowingError

    properties (Access = protected)
        tol = 1e-14;
    end

    properties (Access = private)
        mesh
        symmetricMesh
        scalarField
        symmetricScalarField
    end

    properties (Access = private)
        alpha
        nx
        ny
        x1Max
        x1Min
        x2Max
        x2Min
        symmetricLine
    end

    methods (Access = public)

        function obj = ScalarSymmetrizerTest()
            obj.init();
            obj.createMesh();
            obj.createSymmetricMesh();
            obj.createScalarField();
            obj.computeSymmetricScalarField();
            obj.plotMeshes();
            obj.plotScalarField();
            obj.plotSymmetricScalarField();
        end

    end

    methods (Access = private)

        function init(obj)
            obj.nx = 2;
            obj.ny = 3;
            obj.x1Max = 3;
            obj.x1Min = 1;
            obj.x2Max = 2;
            obj.x2Min = 1;
            obj.alpha = 0;
        end

        function createMesh(obj)
            [x1,x2] = obj.createCoord();
            x3 = zeros(size(x2));
            [F,V] = mesh2tri(x1,x2,x3,'b');
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
            s.symmetricLine = obj.createSymmetricLine();
            mS = Symmetrizer(s);
            m = mS.computeSymmetricMesh();
            obj.symmetricMesh = m;

        end

        function line = createSymmetricLine(obj)
            firstPoint = obj.mesh.coord(1,:);
            line.vector = [cos(obj.alpha);sin(obj.alpha)];
            line.point  = firstPoint;
        end

        function createScalarField(obj)
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            z = x.^2+y.^2;
            obj.scalarField = z;
        end

        function computeSymmetricScalarField(obj)
            s.mesh = obj.mesh;
            s.symmetricLine = obj.createSymmetricLine();
            mS = Symmetrizer(s);
            z  = obj.scalarField;
            zS = mS.symmetrizeScalarField(z);
            obj.symmetricScalarField = zS;
        end

        function plotMeshes(obj)
            figure()
            obj.mesh.plot();
            figure()
            obj.symmetricMesh.plot();
        end

        function plotScalarField(obj)
            z = obj.scalarField();
            m = obj.mesh;
            obj.plotField(m,z);
        end

        function plotSymmetricScalarField(obj)
            z = obj.symmetricScalarField();
            m = obj.symmetricMesh;
            obj.plotField(m,z);
        end

    end

    methods (Access = protected)

        function computeError(obj)
            obj.error = 0;
        end

    end

    methods (Access = private, Static)

        function plotField(mesh,z)
            s.mesh = mesh;
            s.field = z;
            n = NodalFieldPlotter(s);
            n.plot();
        end

    end

end