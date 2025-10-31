classdef GeometricalFunction < handle

    properties (Access = private)
        fHandle
    end

    methods (Access = public)

        function obj = GeometricalFunction(cParams)
            obj.selectHandle(cParams);
        end

        function ls = computeLevelSetFunction(obj,m)
            s.fHandle = obj.fHandle;
            s.ndimf   = 1;
            s.mesh    = m;
            aFun      = AnalyticalFunction(s);
            ls        = aFun.project('P1');
        end

        function f = getHandle(obj)
            f = obj.fHandle;
        end
    end

    methods (Access = private)

        function selectHandle(obj,cParams)
            x1 = @(x) x(1,:,:);
            x2 = @(x) x(2,:,:);
            x3 = @(x) x(3,:,:);
            switch cParams.type
                case 'Square'
                    l  = cParams.length;
                    x0 = cParams.xCoorCenter;
                    y0 = cParams.yCoorCenter;
                    fH = @(x) max(abs(x1(x)-x0),abs(x2(x)-y0))/l - 0.5;
                    obj.fHandle = fH;

                case 'SmoothSquare'
                    l  = cParams.length;
                    x0 = cParams.xCoorCenter;
                    y0 = cParams.yCoorCenter;
                    p  = cParams.pnorm;
                    fH = @(x) (((abs(x1(x)-x0)).^p+(abs(x2(x)-y0)).^p).^(1/p))/l - 0.5;
                    obj.fHandle = fH;

                case 'SquareInclusion'
                    s      = cParams;
                    s.type = 'Square';
                    obj.computeInclusion(s);

                case 'Rectangle'
                   % sx = (cParams.xSide)/2;
                   % sy = (cParams.ySide)/2;
                    % sx = (1-cParams.xSide)/2;
                    % sy = (1-cParams.ySide)/2;    
                    % sx = cos(2*pi*sx);
                    % sy = cos(2*pi*sy);                      
                    % x0 = cParams.xCoorCenter;
                    % y0 = cParams.yCoorCenter;
                    % fH = @(x) max((x1(x)-x0)./(sx),(x2(x)-y0)./(sy)) - 1;
                    % obj.fHandle = fH;

                    sx = cParams.xSide;
                    sy = cParams.ySide;
                    x0 = cParams.xCoorCenter;
                    y0 = cParams.yCoorCenter;
                    fH = @(x) max(abs(x1(x)-x0)./sx,abs(x2(x)-y0)./sy) - 0.5;
                    obj.fHandle = fH;


                case 'RectangleInclusion'
                    s      = cParams;
                    s.type = 'Rectangle';
                    obj.computeInclusion(s);

                case 'SmoothRectangleInclusion'
                    s      = cParams;
                    s.type = 'SmoothRectangle';
                    obj.computeInclusion(s);

                case 'SmoothRectangle'
                     % sx = (cParams.xSide)/2;
                     % sy = (cParams.ySide)/2;
                %   sx = (1-cParams.xSide)/2;
                %   sy = (1-cParams.ySide)/2;    
                %   sx = cos(2*pi*sx);
                %   sy = cos(2*pi*sy);                      
                    % x0 = cParams.xCoorCenter;
                    % y0 = cParams.yCoorCenter;
                    % p  = cParams.pnorm;
                    % fH = @(x) (((x1(x)-x0)./(sx)).^p+((x2(x)-y0)./(sy)).^p).^(1/p) - 1;
                    % obj.fHandle = fH;

                    sx = cParams.xSide;
                    sy = cParams.ySide;
                    x0 = cParams.xCoorCenter;
                    y0 = cParams.yCoorCenter;
                    p  = cParams.pnorm;
                    fH = @(x) ((abs(x1(x)-x0)./sx).^p+(abs(x2(x)-y0)./sy).^p).^(1/p) - 0.5;
                    obj.fHandle = fH;                    

                case 'RectangleRotated'
                    sx = cParams.xSide;
                    sy = cParams.ySide;
                    x0 = cParams.xCoorCenter;
                    y0 = cParams.yCoorCenter;
                    w  = cParams.omegaDeg;
                    R  = [cosd(w),sind(w);-sind(w),cosd(w)];
                    xLoc = @(x) abs(pagemtimes(R,[x1(x)-x0;x2(x)-y0]));
                    fH = @(x) max(pagemtimes([1,0],xLoc(x))/sx,pagemtimes([0,1],xLoc(x))/sy) - 0.5;
                    obj.fHandle = fH;

                case {'Circle','Cylinder'}
                    r  = cParams.radius;
                    x0 = cParams.xCoorCenter;
                    y0 = cParams.yCoorCenter;
                    fH = @(x) (x1(x)-x0).^2+(x2(x)-y0).^2-r^2;
                    obj.fHandle = fH;
                    
                case 'EllipseInclusion'
                    s      = cParams;
                    s.type = 'Ellipse';
                    obj.computeInclusion(s);
                
                case 'Superformula'
                    a  = cParams.semiHorizontalAxis;
                    b  = cParams.semiVerticalAxis;
                    m  = cParams.m;
                    n1 = cParams.n1;
                    n2 = cParams.n2;
                    n3 = cParams.n3;
                    x0 = cParams.xCoorCenter;
                    y0 = cParams.yCoorCenter;
                    
                    phi = @(x) atan2((x2(x)-y0), (x1(x)-x0));
                    r1aux = @(x) abs(cos(m.*phi(x)/4)./a).^n2;
                    r2aux = @(x) abs(sin(m.*phi(x)/4)./b).^n3;
                    r = @(x) (r1aux(x) + r2aux(x)).^(-1./n1);
                    fH = @(x) 1-(((x1(x)-x0)./r(x)).^2 + ((x2(x)-y0)./r(x)).^2); 
                    obj.fHandle = fH;
                    
                case 'SuperformulaInclusion'
                    s      = cParams;
                    s.type = 'Superformula';
                    obj.computeInclusion(s);

                case 'CircleInclusion'
                    s      = cParams;
                    s.type = 'Circle';
                    obj.computeInclusion(s);

                case 'Ellipse'
                    sx = cParams.xSide;
                    sy = cParams.ySide;
                    x0 = cParams.xCoorCenter;
                    y0 = cParams.yCoorCenter;
                    fH = @(x) (((x1(x)-x0).^2)./sx^2)+(((x2(x)-y0).^2)./sy^2) - 1;
                    obj.fHandle = fH;

                case 'Sphere'
                    r  = cParams.radius;
                    x0 = cParams.xCoorCenter;
                    y0 = cParams.yCoorCenter;
                    z0 = cParams.zCoorCenter;
                    fH = @(x) (x1(x)-x0).^2+(x2(x)-y0).^2+(x3(x)-z0).^2-r^2;
                    obj.fHandle = fH;

                case 'SphereInclusion'
                    s      = cParams;
                    s.type = 'Sphere';
                    obj.computeInclusion(s);

                case 'VerticalFiber'
                    l  = cParams.width;
                    x0 = cParams.xCoorCenter;
                    a  = 4/l^2;
                    xm = x0-l/2;
                    xM = x0+l/2;
                    fH = @(x) a*(x1(x)-xm).*(x1(x)-xM);
                    obj.fHandle = fH;

                case 'VerticalNFibers'
                    n    = cParams.nFibers;
                    xmin = cParams.minxCoor;
                    xmax = cParams.maxxCoor;
                    k    = 2*pi*n/(xmax-xmin);
                    fH   = @(x) sin(k*(x1(x)-xmin)+pi/2);
                    obj.fHandle = fH;

                case 'HorizontalFiber'
                    l  = cParams.width;
                    y0 = cParams.yCoorCenter;
                    a  = 4/l^2;
                    ym = y0-l/2;
                    yM = y0+l/2;
                    fH = @(x) a*(x2(x)-ym).*(x2(x)-yM);
                    obj.fHandle = fH;

                case 'HorizontalInclusion'
                    s      = cParams;
                    s.type = 'HorizontalFiber';
                    obj.computeInclusion(s);

                case 'HorizontalNFibers'
                    n    = cParams.nFibers;
                    ymin = cParams.minyCoor;
                    ymax = cParams.maxyCoor;
                    k    = 2*pi*n/(ymax-ymin);
                    fH   = @(x) sin(k*(x2(x)-ymin)+pi/2);
                    obj.fHandle = fH;

                case 'Holes'
                    dim = cParams.dim;
                    n   = cParams.nHoles;
                    l   = cParams.totalLengths;
                    f   = cParams.phases;
                    r   = cParams.phiZero;
                    fH  = @(x) ones(size(x1(x)));
                    for i = 1:dim
                        ni     = n(i);
                        li     = l(i);
                        phasei = f(i);
                        fH     = @(x) fH(x).*cos((ni+1)*x(i,:,:)*pi/li+phasei);
                    end
                    fH = @(x) fH(x)+r-1;
                    obj.fHandle = fH;

                case 'Full'
                    fH = @(x) -1*ones(size(x1(x)));
                    obj.fHandle = fH;

                case 'Given'
                    fH = cParams.fHandle;
                    obj.fHandle = fH;

                case 'Vigdergauz'
                    fH          = LevelSetVigdergauz(cParams);
                    obj.fHandle = fH.getFunctionHandle();

                case 'PeriodicAndOriented'
                    fH          = LevelSetPeriodicAndOriented(cParams);
                    obj.fHandle = fH.getFunctionHandle();
                case 'Hexagon'
                    l  = cParams.radius;
                    n  = cParams.normal;
                    x0 = cParams.xCoorCenter;
                    y0 = cParams.yCoorCenter;
                    p  = 'Inf';
                    fH = @(x) obj.computeHexagonFunction(x,x1,x2,x0,y0,n,p,l);
                    obj.fHandle = fH;
                case 'SmoothHexagon'
                    l  = cParams.radius;
                    n  = cParams.normal;
                    x0 = cParams.xCoorCenter;
                    y0 = cParams.yCoorCenter;
                    p  = cParams.pnorm;
                    fH = @(x) obj.computeHexagonFunction(x,x1,x2,x0,y0,n,p,l);
                    obj.fHandle = fH;
                case 'Circles'
                    r = cParams.r;
                    x0 = cParams.x0;
                    y0 = cParams.y0;
                    fH = @(x) obj.computeCircles(x,x0,y0,r);
                    obj.fHandle = fH;
                case 'FourPerpendicularBars'
                    xL2 = cParams.leftBar_xMax;   % right edge of left bar
                    xR1 = cParams.rightBar_xMin;  % left edge of right bar
                    yB2 = cParams.bottomBar_yMax; % top edge of bottom bar
                    yT1 = cParams.topBar_yMin;    % bottom edge of top bar
                    h = cParams.barWidth;

                    xL1 = xL2 - h;   % left edge of left bar
                    yT2 = yT1 + h;    % top edge of top bar
                    yB1 = yB2 - h; % bottom edge of bottom bar
                    xR2 = xR1 + h;  % right edge of right bar

                    fV1 = @(x) max( xL1 - x1(x), x1(x) - xL2 );
                    fV2 = @(x) max( xR1 - x1(x), x1(x) - xR2 );
                    fH1 = @(x) max( yB1 - x2(x), x2(x) - yB2 );
                    fH2 = @(x) max( yT1 - x2(x), x2(x) - yT2 );

                    fH = @(x) min( min(fV1(x), fV2(x)), min(fH1(x), fH2(x)) );

                    obj.fHandle = fH;


            end
        end

        function computeInclusion(obj,s)
            obj.selectHandle(s);
            fH          = obj.fHandle;
            obj.fHandle = @(x) -fH(x);
        end

    end

    methods (Access = private, Static)

        function d = computeHexagonFunction(x,x1,x2,x0,y0,n,p,l)
            vx     = x1(x)-x0;
            vy     = x2(x)-y0;
            nS     = size(n,1);
            nGauss = size(x,2);
            nElem  = size(x,3);
            vn = zeros(1,nGauss,nElem,nS);
            for i = 1:nS
                nx = n(i,1);
                ny = n(i,2);
                vn(:,:,:,i) = abs(vx*nx + vy*ny);
            end
            normVn = vecnorm(vn,p,4);
            d = (normVn/(l*(sqrt(3)/2)))-1;
        end

        function fH = computeCircles(x,x0,y0,r)
            n = length(r);
            for i =1:n
                f(:,:,:,i) = (x(1,:,:)-x0(i)).^2+(x(2,:,:)-y0(i)).^2-r(i)^2;
            end   
            fH = min(f,[],4);
        end

    end
end