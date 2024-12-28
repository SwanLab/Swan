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
            obj.createOrientations();
            obj.dehomogenize();
        end

    end

    methods (Access = private)

        function init(obj)
            obj.meshSize = 0.03;
            obj.nCells   = [10 10; 20 20];
            obj.xmin = 0;
            obj.xmax = 2;
            obj.ymin = 0;
            obj.ymax = 1;
            %obj.widthH = 0.4;
            obj.widthH = 0.85;
            obj.widthW = 0.85;
        end

        function createMesh(obj)
            h = obj.meshSize;
            xv = obj.xmin:h:obj.xmax;
            yv = obj.ymin:h:obj.ymax;
            [X,Y] = meshgrid(xv,yv);
            s.coord(:,1) = X(:);
            s.coord(:,2) = Y(:);
            [F,V] = mesh2tri(X,Y,zeros(size(X)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;

            m = Mesh.create(s);
            obj.mesh = m;
        end

        function createOrientations(obj)
            coord = obj.mesh.coord;
            x1 = coord(:,1);
            x2 = coord(:,2);
            x10 = (max(x1(:))+min(x1(:)))/2;
            x20 = -0.1*max(x2(:));                                    
            r = sqrt((x1-x10).^2+(x2-x20).^2);
            fR = obj.normalize([(x1-x10)./r,(x2-x20)./r]);            
            fT = obj.normalize([-(x2-x20)./r,(x1-x10)./r]);
            obj.orientation{1} = obj.createOrientationField(fR);
            obj.orientation{2} = obj.createOrientationField(fT);
        end 

        function fN = normalize(obj,f)
            fN = f./vecnorm(f,2,2);
        end

        function f = createOrientationField(obj,fV)
            fD = obj.createP1DiscontinousOrientation(fV);
            f  = obj.computeAngleOrientation(fD);            
        end

     function fS = computeAngleOrientation(obj,f)
            aV = f.fValues;
            aV1 = aV(:,1);
            aV2 = aV(:,2);
            fN(:,1) = aV1.^2-aV2.^2;
            fN(:,2) = 2*aV1.*aV2;
            s.fValues = fN;
            s.mesh    = obj.mesh;
            s.order   = 'P1D';
            s.ndimf   = 2;
            fS = LagrangianFunction(s);            
        end

        function aF = createP1DiscontinousOrientation(obj,fV)
            s.fValues = fV;
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            s.ndimf   = 2;
            aF = LagrangianFunction(s);
            aF = project(aF,'P1D');
        end

        function plotOrientation(obj)            
            plotVector(obj.orientation{1},4);
            plotVector(obj.orientation{2},4);
        end

        function s = createLevelSetCellParams(obj)
            s.type  = 'RectangleInclusion';
            s.xSide = obj.createFunction(obj.widthH,1);
            s.ySide = obj.createFunction(obj.widthW,2);
            s.ndim   = 2;
        end

        function f = createFunction(obj,value,dir)
            s.fHandle = @(x) obj.variationFunction(value,x,dir);
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            f = AnalyticalFunction(s);
            f = f.project('P1D');
        end

        function f = variationFunction(obj,mV,x,dir)
            xV = x(dir,:,:);
            I  = ones(size(xV));
            xmin = min(xV(:));
            xmax = max(xV(:));
            incX   = (xmax-xmin);
            xM     = (xmax+xmin)/2;
            mMin = 0.01;
            mMax = 0.99;
            incM = mMax - mMin;
            tMax_min = (mV - mMin)/(2*incM)*incX/(xM - xmin);
            tMax_max = (mMax - mV)/(2*incM)*incX/(xmax - xM);

            %t = min(tMax_min,tMax_max);

            t = 0;

            f = mV + 2*t*(incM)/(incX)*(xV-xM);
        end

        function dehomogenize(obj)
            s.nCells             = obj.nCells;
            s.cellLevelSetParams = obj.createLevelSetCellParams();
            s.mesh               = obj.mesh;
            s.orientation        = obj.orientation;
            d = Dehomogenizer(s);
            ls = d.compute();
            d.plot();
            obj.levelSet = ls;
        end

    end

end
