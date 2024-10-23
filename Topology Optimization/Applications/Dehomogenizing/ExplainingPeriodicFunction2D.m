classdef ExplainingPeriodicFunction2D < handle

    properties (Access = private)
        mesh
        orientation
        levelSet
    end

    properties (Access = private)
        meshSize
        xmin
        xmax
        ymin
        ymax
        widthH
        widthW
        nCells
        alpha
    end

    methods (Access = public)

        function obj = ExplainingPeriodicFunction2D()
            obj.init();
            obj.createMesh();
            obj.createAngle();
            obj.createOrientation();
            obj.dehomogenize();
        end

    end

    methods (Access = private)

        function init(obj)
            obj.meshSize = 0.01;
            obj.nCells   = 5;
            obj.xmin = 0;
            obj.xmax = 1;
            obj.ymin = 0;
            obj.ymax = 1;
            obj.widthH = 0.9;
            obj.widthW = 0.9;
        end

        function createMesh(obj)
            h = obj.meshSize;
            xv = obj.xmin:h:obj.xmax;
            yv = obj.ymin:h:obj.ymax;
            [X,Y] = meshgrid(xv,yv);
            s.coord(:,1) = X(:);
            s.coord(:,2) = Y(:);
            [F,V] = mesh2tri(X,Y,zeros(size(X)),'f');
            s.coord  = V(:,1:2);
            s.connec = F;

            m = Mesh.create(s);
            obj.mesh = m;
        end

        function createAngle(obj)
            m = obj.mesh;
            x1 = m.coord(:,1);
            beta = zeros(size(x1));
            alpha = beta/2;
            s.fValues = alpha;
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            aF = LagrangianFunction(s);
            obj.alpha = aF;
        end


        function createOrientation(obj)
            al = obj.alpha.fValues;
            a = [cos(al), sin(al)];
            s.fValues = a;
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            aF = LagrangianFunction(s);            
            obj.orientation{1} = aF;
            a = [-sin(al), cos(al)];
            s.fValues = a;
            s.mesh    = obj.mesh;
            aF = LagrangianFunction(s);            
            obj.orientation{2} = aF;
        end

        function plotOrientation(obj)
            figure()
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            t  = obj.orientation{1}.fValues;
            ct = (t(:,1));
            st = (t(:,2));
            quiver(x,y,ct,st)
        end

        function s = createLevelSetCellParams(obj)            
            s.type  = 'RectangleInclusion';
            s.xSide = obj.createConstantFunction(obj.widthH);
            s.ySide = obj.createConstantFunction(obj.widthW);
            s.ndim   = 2;
        end

        function f = createConstantFunction(obj,value)
            s.fHandle = @(x) value*ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            f = AnalyticalFunction(s);
            f = f.project('P1D');
        end

        function dehomogenize(obj)
            s.nCells             = obj.nCells;
            s.cellLevelSetParams = obj.createLevelSetCellParams();
            s.mesh               = obj.mesh;
            s.theta              = obj.orientation;%.project('P0');%atan2(obj.orientation(:,2),obj.orientation(:,1));
            d = Dehomogenizer(s);
            ls = d.compute();
            d.plot();
            obj.levelSet = ls;
        end

    end

end
