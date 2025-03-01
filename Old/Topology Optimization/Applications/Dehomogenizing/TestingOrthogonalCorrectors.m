classdef TestingOrthogonalCorrectors < handle

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (Access = private)
        mesh
        orientation
        interpolator
        orthogonalCorrector
    end

    methods (Access = public)

        function obj = TestingOrthogonalCorrectors()
            %  obj.createToyMesh();
            %  obj.createToyOrientation();
            obj.createBenchmarkMesh();
            obj.createBenchmarkOrientation();
            obj.plotOrientationVector();
            obj.createInterpolator();
            obj.createOrthogonalCorrector();
        end

    end

    methods (Access = private)

        function createToyMesh(obj)
            x = [66.90 128.89 115.25 76.73 26.84 157.24 ...
                168.58 141.74 98.65 45.74 3.78 2.65 31.37 ...
                83.53 137.21 174.63 147.03 99.03 38.55];
            y = [89.20 89.58 120.20 130.40 119.44 113.39 ...
                153.08 154.59 160.26 164.80 152.32 220.74 ...
                207.13 202.977 190.88 192.77 239.64 241.90 242.28];
            s.coord(:,1) = x;
            s.coord(:,2) = y;
            s.connec = delaunay(s.coord);
            m = Mesh.create(s);
            m.plot();
            obj.mesh = m;
        end

        function createToyOrientation(obj)
            alpha = pi/180*[190 140 300 150 185 ...
                120 75  45  150 130 ...
                160+180 0   190  5  -5   ...
                170 0   190 190];
            a(:,1) = cos(alpha);
            a(:,2) = sin(alpha);
            obj.orientation = a;
        end

        function createBenchmarkMesh(obj)
            h = 0.03;
            xmin = 0.50;
            xmax = 2.0;
            ymin = 0.25;
            ymax = 1.75;
            xv = xmin:h:xmax;
            yv = ymin:h:ymax;
            [X,Y] = meshgrid(xv,yv);
            s.coord(:,1) = X(:);
            s.coord(:,2) = Y(:);
            s.connec = delaunay(s.coord);
            m = Mesh.create(s);
            obj.mesh = m;
        end

        function createBenchmarkOrientation(obj)
            s1 = 0.32;
            s2 = -0.8;
            x1 = obj.mesh.coord(:,1);
            x2 = obj.mesh.coord(:,2);
            v(:,1) = cos(pi*(x1 + s1*x2));
            v(:,2) = cos(pi*(x2 + s2*x1));
            beta = atan2(v(:,2),v(:,1));
            alpha = beta/2;
            obj.orientation(:,1) = cos(alpha);
            obj.orientation(:,2) = sin(alpha);
        end

        function plotOrientationVector(obj)
            %figure()
            a = obj.orientation;
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            ax = a(:,1);
            ay = a(:,2);
            quiver(x,y,ax,ay);
        end

        function createInterpolator(obj)
            s.meshCont    = obj.mesh;
            s.meshDisc    = obj.mesh.createDiscontinousMesh;
            s.orientation = obj.orientation;
            s = SymmetricContMapCondition(s);
            sC = s.computeCondition();
            obj.interpolator = sC;
        end

        function createOrthogonalCorrector(obj)
            s.mesh            = obj.mesh;
            s.interpolator    = obj.interpolator;
            s.correctorValue  = obj.computeCorrector();
            o = OrthogonalCorrectorComputer(s);
            oC = o.compute();
            o.plot();
            obj.orthogonalCorrector = oC;
        end

        function cV = computeCorrector(obj)
            s.mesh               = obj.mesh;
            s.orientation        = obj.orientation;
            s.singularityCoord   = obj.computeSingularities();
            c = CorrectorComputer(s);
            cV = c.compute();
            c.plot()
        end

        function sCoord = computeSingularities(obj)
            s.mesh        = obj.mesh;
            s.orientation = obj.orientation;
            sF = SingularitiesFinder(s);
            isS = sF.computeSingularElements();
            sF.plot();
            coordB = obj.mesh.computeBaricenter();
            coordB = transpose(coordB);
            sCoord =  coordB(isS,:);
            sCoord = sCoord(1,:);
        end

    end



end