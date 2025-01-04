classdef ExplainingPeriodicFunction2D < handle

    properties (Access = private)
        mesh
        orientationA
        orientationB
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
    end

    methods (Access = public)

        function obj = ExplainingPeriodicFunction2D()
            obj.init();
            obj.createMesh();
            obj.createOrientationA();
            obj.createOrientationB();
            obj.dehomogenize();
        end

    end

    methods (Access = private)

        function init(obj)
            obj.meshSize = 0.03;
            %obj.meshSize = 0.01;
            n0 = 4; nF = 30; nN = 3;
            obj.nCells   = round(repmat(linspace(n0,nF,nN),2,1)'/2)*2;
             obj.xmin = 0;
             obj.xmax = 2;
             obj.ymin = 0;
             obj.ymax = 1;
           %obj.xmin = -0.005;
           %obj.xmax = 0.005;
           %obj.ymin = 0;
           %obj.ymax = 0.01;
            
            obj.widthH = 0.85;%0.85;
            obj.widthW = 0.85;%0.85;
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

        function createOrientationA(obj)
            coord = obj.mesh.coord;
            x1 = coord(:,1);
            x2 = coord(:,2);
            x10 = (max(x1(:))+min(x1(:)))/2+0;
            x20 = -0.1*max(x2(:));                                    
            r = sqrt((x1-x10).^2+(x2-x20).^2);
            fR = obj.normalize([(x1-x10)./r,(x2-x20)./r]);            
            fT = obj.normalize([-(x2-x20)./r,(x1-x10)./r]);
            obj.orientationA{1} = obj.createP1DiscontinousOrientation(fR);
            obj.orientationA{2} = obj.createP1DiscontinousOrientation(fT);
        end 

        function fN = normalize(obj,f)
            fN = f./vecnorm(f,2,2);
        end

        function aF = createP1DiscontinousOrientation(obj,fV)
            s.fValues = fV;
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            s.ndimf   = 2;
            aF = LagrangianFunction(s);
            aF = project(aF,'P1D');
        end

      function createOrientationB(obj)
            for iDim = 1:obj.mesh.ndim
                aI = obj.orientationA{iDim};
                bI = obj.computeDoubleOrientationAngle(aI);
                obj.orientationB{iDim} = bI;
            end
        end

     function fS = computeDoubleOrientationAngle(obj,f)
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

     function fS = createP1DiscontinousFunction(obj,f)
            s.fValues = f;
            s.mesh    = obj.mesh;
            s.order   = 'P1D';
            s.ndimf   = 2;
            fS = LagrangianFunction(s);            
        end       

        function plotOrientation(obj)            
            plotVector(obj.orientationA{1},4);
            plotVector(obj.orientationA{2},4);
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

        function a = createOrientationAfromB(obj)
            b = obj.orientationB{1};
            beta = atan2(b.fValues(:,2),b.fValues(:,1));
          %  alpha = beta/2;%0*ones(size(beta));%beta/2;            
            alpha = 0*ones(size(beta));%beta/2;            
            a1 = [cos(alpha), sin(alpha)];
            a2 = [-sin(alpha), cos(alpha)];

           % a1([1 3 4],1) = 0;
           % a1([1 3 4],2) = -1;
           % a1([2 5 6],1) = 0;
           % a1([2 5 6],2) = 1;
           % a2([1 3 4],1) = 1;
           % a2([1 3 4],2) = 0;
           % a2([2 5 6],1) = -1;
           % a2([2 5 6],2) = 0;                       


            a{1} = obj.createP1DiscontinousFunction(a1);
            a{2} = obj.createP1DiscontinousFunction(a2);
        end


        function dehomogenize(obj)
            s.nCells             = obj.nCells;
            s.cellLevelSetParams = obj.createLevelSetCellParams();
            s.mesh               = obj.mesh;
        %  s.orientationA      = obj.orientationA;
            s.orientationA       = obj.createOrientationAfromB(); 
            d = Dehomogenizer(s);
            ls = d.compute();
            d.plot();
            obj.levelSet = ls;
        end

    end

end
