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
                    fH = @(x) ((abs(x1(x)-x0).^p+abs(x2(x)-y0).^p).^(1/p))/l - 0.5;
                    obj.fHandle = fH;

                case 'SquareInclusion'
                    s      = cParams;
                    s.type = 'Square';
                    obj.computeInclusion(s);

                case 'Rectangle'
                    sx = cParams.xSide;
                    sy = cParams.ySide;
                    x0 = cParams.xCoorCenter;
                    y0 = cParams.yCoorCenter;
                    fH = @(x) max(abs(x1(x)-x0)/sx,abs(x2(x)-y0)/sy) - 0.5;
                    obj.fHandle = fH;

                case 'RectangleInclusion'
                    s      = cParams;
                    s.type = 'Rectangle';
                    obj.computeInclusion(s);

                case 'SmoothRectangle'
                    sx = cParams.xSide;
                    sy = cParams.ySide;
                    x0 = cParams.xCoorCenter;
                    y0 = cParams.yCoorCenter;
                    p  = cParams.pnorm;
                    fH = @(x) ((abs(x1(x)-x0)/sx).^p+(abs(x2(x)-y0)/sy).^p).^(1/p) - 0.5;
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

                case 'CircleInclusion'
                    s      = cParams;
                    s.type = 'Circle';
                    obj.computeInclusion(s);

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
                    vig         = LevelSetVigdergauz(cParams);
                    obj.fHandle = vig.getFunctionHandle();

                case 'PeriodicAndOriented'
                    perOr       = LevelSetPeriodicAndOriented(cParams);
                    obj.fHandle = perOr.getFunctionHandle();
                case 'Custom'
                    r    = 0.1;
                    fH11 = @(x) (x1(x)-0.25).^2+(x2(x)-0.25).^2+(x3(x)-0.25).^2-r^2 <= 0;
                    fH12 = @(x) (x1(x)-0.25).^2+(x2(x)-0.25).^2+(x3(x)-0.50).^2-r^2 <= 0;
                    fH13 = @(x) (x1(x)-0.25).^2+(x2(x)-0.25).^2+(x3(x)-0.75).^2-r^2 <= 0;
                    fH14 = @(x) (x1(x)-0.25).^2+(x2(x)-0.50).^2+(x3(x)-0.25).^2-r^2 <= 0;
                    fH15 = @(x) (x1(x)-0.25).^2+(x2(x)-0.50).^2+(x3(x)-0.50).^2-r^2 <= 0;
                    fH16 = @(x) (x1(x)-0.25).^2+(x2(x)-0.50).^2+(x3(x)-0.75).^2-r^2 <= 0;
                    fH17 = @(x) (x1(x)-0.25).^2+(x2(x)-0.75).^2+(x3(x)-0.25).^2-r^2 <= 0;
                    fH18 = @(x) (x1(x)-0.25).^2+(x2(x)-0.75).^2+(x3(x)-0.50).^2-r^2 <= 0;
                    fH19 = @(x) (x1(x)-0.25).^2+(x2(x)-0.75).^2+(x3(x)-0.75).^2-r^2 <= 0;
                    fH21 = @(x) (x1(x)-0.50).^2+(x2(x)-0.25).^2+(x3(x)-0.25).^2-r^2 <= 0;
                    fH22 = @(x) (x1(x)-0.50).^2+(x2(x)-0.25).^2+(x3(x)-0.50).^2-r^2 <= 0;
                    fH23 = @(x) (x1(x)-0.50).^2+(x2(x)-0.25).^2+(x3(x)-0.75).^2-r^2 <= 0;
                    fH24 = @(x) (x1(x)-0.50).^2+(x2(x)-0.50).^2+(x3(x)-0.25).^2-r^2 <= 0;
                    fH25 = @(x) (x1(x)-0.50).^2+(x2(x)-0.50).^2+(x3(x)-0.50).^2-r^2 <= 0;
                    fH26 = @(x) (x1(x)-0.50).^2+(x2(x)-0.50).^2+(x3(x)-0.75).^2-r^2 <= 0;
                    fH27 = @(x) (x1(x)-0.50).^2+(x2(x)-0.75).^2+(x3(x)-0.25).^2-r^2 <= 0;
                    fH28 = @(x) (x1(x)-0.50).^2+(x2(x)-0.75).^2+(x3(x)-0.50).^2-r^2 <= 0;
                    fH29 = @(x) (x1(x)-0.50).^2+(x2(x)-0.75).^2+(x3(x)-0.75).^2-r^2 <= 0;
                    fH31 = @(x) (x1(x)-0.75).^2+(x2(x)-0.25).^2+(x3(x)-0.25).^2-r^2 <= 0;
                    fH32 = @(x) (x1(x)-0.75).^2+(x2(x)-0.25).^2+(x3(x)-0.50).^2-r^2 <= 0;
                    fH33 = @(x) (x1(x)-0.75).^2+(x2(x)-0.25).^2+(x3(x)-0.75).^2-r^2 <= 0;
                    fH34 = @(x) (x1(x)-0.75).^2+(x2(x)-0.50).^2+(x3(x)-0.25).^2-r^2 <= 0;
                    fH35 = @(x) (x1(x)-0.75).^2+(x2(x)-0.50).^2+(x3(x)-0.50).^2-r^2 <= 0;
                    fH36 = @(x) (x1(x)-0.75).^2+(x2(x)-0.50).^2+(x3(x)-0.75).^2-r^2 <= 0;
                    fH37 = @(x) (x1(x)-0.75).^2+(x2(x)-0.75).^2+(x3(x)-0.25).^2-r^2 <= 0;
                    fH38 = @(x) (x1(x)-0.75).^2+(x2(x)-0.75).^2+(x3(x)-0.50).^2-r^2 <= 0;
                    fH39 = @(x) (x1(x)-0.75).^2+(x2(x)-0.75).^2+(x3(x)-0.75).^2-r^2 <= 0;
                    in   = @(x) fH11(x) | fH12(x) | fH13(x) | fH14(x) | fH15(x) | fH16(x) | fH17(x) | fH18(x) | fH19(x) ...
                                 | fH21(x) | fH22(x) | fH23(x) | fH24(x) | fH25(x) | fH26(x) | fH27(x) | fH28(x) | fH29(x) ...
                                  | fH31(x) | fH32(x) | fH33(x) | fH34(x) | fH35(x) | fH36(x) | fH37(x) | fH38(x) | fH39(x);
                    fH   = @(x) in(x).*1 + not(in(x)).*-1;
                    obj.fHandle = fH;
            end
        end

        function computeInclusion(obj,s)
            obj.selectHandle(s);
            fH          = obj.fHandle;
            obj.fHandle = @(x) -fH(x);
        end
    end
end