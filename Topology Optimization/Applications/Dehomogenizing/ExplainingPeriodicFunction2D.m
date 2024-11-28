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
            obj.meshSize = 0.03;
            obj.nCells   = [10 10; 20 20];
            obj.xmin = 0;
            obj.xmax = 2;
            obj.ymin = 0;
            obj.ymax = 1;
            %obj.widthH = 0.4;
            obj.widthH = 0.8;
            obj.widthW = 0.8;
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
            s.fHandle = @(x) obj.createAlphaValues(x);
            s.mesh    = obj.mesh;
            s.ndimf   = 1;
            obj.alpha = AnalyticalFunction(s);


            gradA = Grad(obj.alpha.project('P1'));
            gradA = gradA.project('P1',obj.mesh);
            t = Divergence(gradA);
            t.project('P1',obj.mesh).plot()
        end

        function f = createAlphaValues(obj,x)
            x1 = x(1,:,:);
            x2 = x(2,:,:);
            x10 = (max(x1(:))+min(x1(:)))/2;
            x20 = 0;            
            f = atan2(x2-x20 +0.5*(max(x2(:))),x1-x10);
            isLeft = x1 < (min(x1(:))+ max(x1(:)))/2;
            f(isLeft) = f(isLeft) + pi;
        end

        function createOrientation(obj)
            nDim = obj.mesh.ndim;
            for iDim = 1:nDim
                s.operation = @(xV) obj.createOrientationFunction(iDim,xV);
                s.ndimf     = 2;
                aF = DomainFunction(s);
                aF = aF.project('P1',obj.mesh);
                obj.orientation{iDim} = aF;
            end

        end

        function or = createOrientationFunction(obj,iDim,xV)
            alphaV = obj.alpha.evaluate(xV);
            if iDim == 1                
                or(1,:,:) = cos(alphaV);
                or(2,:,:) = sin(alphaV);
            else
                alphaV = obj.alpha.evaluate(xV);
                or(1,:,:) = -sin(alphaV);
                or(2,:,:) = cos(alphaV);
            end
        end

        function plotOrientation(obj,varargin)
            obj.orientation{1}.project('P1',obj.mesh).plotVector();
            obj.orientation{2}.project('P1',obj.mesh).plotVector();
        end

        function s = createLevelSetCellParams(obj)
            s.type  = 'RectangleInclusion';
            s.xSide = obj.createFunction(obj.widthH,1);
            s.ySide = obj.createFunction(obj.widthW,2);
            s.ndim   = 2;
        end

        function f = createFunction(obj,value,dir)
            s.fHandle = @(x) obj.variationFunction(value,x,dir);
            %  s.fHandle = @(x) value*ones(size(x(1,:,:)));%x(1,:,:);%ones(size(x(1,:,:)));
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
            s.theta              = obj.orientation;%.project('P0');%atan2(obj.orientation(:,2),obj.orientation(:,1));
            d = Dehomogenizer(s);
            ls = d.compute();
            d.plot();
            obj.levelSet = ls;
        end

    end

end
