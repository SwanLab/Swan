classdef StructuralComputer < handle

    properties (Access = private)
        mesh
        orientation
        levelSet
    end

    properties (Access = private)
        testName
        meshSize
        singularitiesData
        xmin
        xmax
        ymin
        ymax
        widthH
        widthW
        nCells
    end

    methods (Access = public)

        function obj = StructuralComputer()
            obj.init();
            obj.createMesh();
            obj.createOrientation();
            obj.dehomogenize();
        end

    end

    methods (Access = private)

        function init(obj)
            obj.meshSize = 0.05;%0.0221;%0.09;%0.0221;%0.0521 %0.0221;0.0921
            obj.nCells   = [20 20];%linspace(60,62,40);%45;   %45
            obj.xmin = 0.5;
            obj.xmax = 2.0;
            obj.ymin = 0.25;
            obj.ymax = 1.75;
            obj.singularitiesData = [0.32,-0.8];
            obj.widthH = 0.87;
            obj.widthW = 0.87;
        end

        function createMesh(obj)
            h = obj.meshSize;
            xv = obj.xmin:h:obj.xmax;
            yv = obj.ymin:h:obj.ymax;
            [X,Y] = meshgrid(xv,yv);
            s.coord(:,1) = X(:);
            s.coord(:,2) = Y(:);
            s.connec = delaunay(s.coord);
            m = Mesh(s);
            obj.mesh = m;
        end


        function createOrientation(obj)
            m = obj.mesh;
            x1 = m.coord(:,1);
            x2 = m.coord(:,2);
            theta = 0*atan(x1./x2);
            obj.orientation(:,1) = cos(theta);
            obj.orientation(:,2) = sin(theta);
        end

        function plotOrientation(obj)
            figure()
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            t  = obj.orientation;
            ct = cos(t(:,1));
            st = sin(t(:,1));
            quiver(x,y,ct,st)
        end

        function s = createLevelSetCellParams(obj)
            I        = ones(size(obj.mesh.coord,1),1);
            s.type   = 'rectangleInclusion';
            s.widthH = obj.widthH*I;
            s.widthV = obj.widthW*I;
            s.ndim   = 2;
        end

        function dehomogenize(obj)
            s.nCells             = obj.nCells;
            s.cellLevelSetParams = obj.createLevelSetCellParams();
            s.mesh               = obj.mesh;
            s.theta              = atan2(obj.orientation(:,2),obj.orientation(:,1));
            d = Dehomogenizer(s);
            ls = d.compute();
            d.plot();
            obj.levelSet = ls;
        end

    end

end