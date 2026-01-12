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

                case 'ThreeRectangles'
                    sx1 = cParams.xSide1;
                    sy1 = cParams.ySide1;
                    x01 = cParams.xCoorCenter1;
                    y01 = cParams.yCoorCenter1;
                    sx2 = cParams.xSide2;
                    sy2 = cParams.ySide2;
                    x02 = cParams.xCoorCenter2;
                    y02 = cParams.yCoorCenter2;
                    sx3 = cParams.xSide3;
                    sy3 = cParams.ySide3;
                    x03 = cParams.xCoorCenter3;
                    y03 = cParams.yCoorCenter3;
                    fH = @(x) min(min(max(abs(x1(x)-x01)./sx1,abs(x2(x)-y01)./sy1), max(abs(x1(x)-x02)./sx2,abs(x2(x)-y02)./sy2)), max(abs(x1(x)-x03)./sx3,abs(x2(x)-y03)./sy3)) - 0.5;
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
                                        

                case 'FourPerpendicularBarsWithCrack'
                    xL2 = cParams.leftBar_xMax;   % right edge of left bar
                    xR1 = cParams.rightBar_xMin;  % left edge of right bar
                    yB2 = cParams.bottomBar_yMax; % top edge of bottom bar
                    yT1 = cParams.topBar_yMin;    % bottom edge of top bar
                    h = cParams.barWidth;
                    hCrack = cParams.hCrack;   % vertical thickness of the crack

                    xL1 = xL2 - h;   % left edge of left bar
                    yT2 = yT1 + h;    % top edge of top bar
                    yB1 = yB2 - h; % bottom edge of bottom bar
                    xR2 = xR1 + h;  % right edge of right bar

                    fV1 = @(x) max( xL1 - x1(x), x1(x) - xL2 );
                    fV2 = @(x) max( xR1 - x1(x), x1(x) - xR2 );
                    fH1 = @(x) max( yB1 - x2(x), x2(x) - yB2 );
                    fH2 = @(x) max( yT1 - x2(x), x2(x) - yT2 );
                
                    fBars = @(x) min( min(fV1(x), fV2(x)), min(fH1(x), fH2(x)) );
                
                    y0_crack = (yB2 + yT1) / 2;   % vertical center of crack
                    yC1 = y0_crack - hCrack/2;
                    yC2 = y0_crack + hCrack/2;
                    delta = 0.05 * (xR2 - xR1);
                    fCrack = @(x) max( max(xR1 - delta - x1(x), x1(x) - xR2 - delta), ...
                                       max(yC1 - x2(x), x2(x) - yC2) );
                
                    fH = @(x) max( fBars(x), -fCrack(x) );
                    obj.fHandle = fH;

                case 'Prism'
                    sx = cParams.xSide;
                    sy = cParams.ySide;
                    sz = cParams.zSide;
                    x0 = cParams.xCoorCenter;
                    y0 = cParams.yCoorCenter;
                    z0 = cParams.zCoorCenter;
                    fH = @(x) max(max(abs(x1(x)-x0)./sx,abs(x2(x)-y0)./sy), abs(x3(x)-z0)./sz) - 0.5;
                    obj.fHandle = fH;


                case 'ThreePrisms'
                    sx1 = cParams.xSide1;
                    sy1 = cParams.ySide1;
                    sz1 = cParams.zSide1;
                    x01 = cParams.xCoorCenter1;
                    y01 = cParams.yCoorCenter1;
                    z01 = cParams.zCoorCenter1;
                    sx2 = cParams.xSide2;
                    sy2 = cParams.ySide2;
                    sz2 = cParams.zSide2;
                    x02 = cParams.xCoorCenter2;
                    y02 = cParams.yCoorCenter2;
                    z02 = cParams.zCoorCenter2;     
                    sx3 = cParams.xSide3;
                    sy3 = cParams.ySide3;
                    sz3 = cParams.zSide3;
                    x03 = cParams.xCoorCenter3;
                    y03 = cParams.yCoorCenter3;
                    z03 = cParams.zCoorCenter3;
                    fH = @(x) min(min(max(max(abs(x1(x)-x01)./sx1,abs(x2(x)-y01)./sy1), abs(x3(x)-z01)./sz1), max(max(abs(x1(x)-x02)./sx2,abs(x2(x)-y02)./sy2), abs(x3(x)-z02)./sz2)), max(max(abs(x1(x)-x03)./sx3,abs(x2(x)-y03)./sy3), abs(x3(x)-z03)./sz3)) - 0.5;
                    obj.fHandle = fH;    

                case 'ThreeRectanglesInclusion'
                    s      = cParams;
                    s.type = 'ThreeRectangles';
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

                case 'TwoCircles'
                    r  = cParams.radius;
                    x0 = cParams.xCoorCenter;
                    y0 = cParams.yCoorCenter;
                    x02 = cParams.xCoorCenter2;
                    y02 = cParams.yCoorCenter2;
                    fH = @(x) min((x1(x)-x0).^2+(x2(x)-y0).^2-r^2, (x1(x)-x02).^2+(x2(x)-y02).^2-r^2);
                    obj.fHandle = fH;

                                     
                case 'TwoCirclesInclusion'
                    s      = cParams;
                    s.type = 'TwoCircles';
                    obj.computeInclusion(s);


                case 'RingSDF'
                    Rin = cParams.innerRadius;
                    Rout = cParams.outerRadius;
                    x0  = cParams.xCoorCenter;
                    y0  = cParams.yCoorCenter;
                
                    if Rout <= Rin
                        error('outerRadius must be greater than innerRadius');
                    end
                
                    fH = @(x) max( Rin - sqrt( (x1(x)-x0).^2 + (x2(x)-y0).^2 ), ...
                                    sqrt( (x1(x)-x0).^2 + (x2(x)-y0).^2 ) - Rout );
                
                    % Inside the ring: fH(x) < 0; magnitude = distance to nearest boundary
                    obj.fHandle = fH;

                case 'RingInclusion'
                    s      = cParams;
                    s.type = 'RingSDF';
                    obj.computeInclusion(s);

                case 'RingWithHorizontalCrack'
                    Rin = cParams.innerRadius;
                    Rout = cParams.outerRadius;
                    x0  = cParams.xCoorCenter;
                    y0  = cParams.yCoorCenter;
                    w = cParams.crackWidth;
                
                    if Rout <= Rin
                        error('outerRadius must be greater than innerRadius');
                    end
                
                    % Helpers
                    r  = @(x) sqrt( (x1(x)-x0).^2 + (x2(x)-y0).^2 );
                
                    % Annulus SDF: negative inside the ring
                    phi_ring  = @(x) max(Rin - r(x), r(x) - Rout);
                
                    % Crack SDF (intersection of: horizontal strip, half-plane x>=x0, and ring radii)
                    phi_crack = @(x) max( max(abs(x2(x)-y0) - w/2, -(x1(x)-x0)), ...
                                          max(r(x) - Rout, Rin - r(x)) );
                
                    % Set difference: Ring \ Crack  =>  max(φ_ring, -φ_crack)
                    fH = @(x) max( phi_ring(x), -phi_crack(x) );
                
                    obj.fHandle = fH;

                case 'RingWithHorizontalCrackInclusion'
                    s      = cParams;
                    s.type = 'RingWithHorizontalCrack';
                    obj.computeInclusion(s);

                case 'FiveCircles'
                    r  = cParams.radius;
                    Lx = cParams.width;         % total width of domain
                    Ly = cParams.height;        % total height of domain
                
                    % Corner and center circle coordinates
                    centers = [ ...
                        0,   0;    % bottom-left
                        Lx,  0;    % bottom-right
                        0,   Ly;   % top-left
                        Lx,  Ly;   % top-right
                        Lx/2, Ly/2 % center
                    ];
                
                    fH = @(x) inf;
                    for i = 1:size(centers,1)
                        xc = centers(i,1);
                        yc = centers(i,2);
                        fH = @(x) min(fH(x), (x1(x)-xc).^2 + (x2(x)-yc).^2 - r^2);
                    end
                    obj.fHandle = fH;

                case 'FiveCirclesInclusion'
                    s      = cParams;
                    s.type = 'FiveCircles';
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

                case 'TwoSpheres'
                    r  = cParams.radius;
                    x0 = cParams.xCoorCenter;
                    y0 = cParams.yCoorCenter;
                    z0 = cParams.zCoorCenter;
                    x02 = cParams.xCoorCenter2;
                    y02 = cParams.yCoorCenter2;
                    z02 = cParams.zCoorCenter2;
                    fH = @(x) min((x1(x)-x0).^2+(x2(x)-y0).^2+(x3(x)-z0).^2-r^2, (x1(x)-x02).^2+(x2(x)-y02).^2+(x3(x)-z02).^2-r^2);
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

                case 'PerperndicularNFiber'
                    n    = cParams.nFibers;
                    xmin = cParams.minxCoor;
                    xmax = cParams.maxxCoor;
                    k    = 2*pi*n/(xmax-xmin);
                    fV   = @(x) sin(k*(x1(x)-xmin)+pi/2);

                    n    = cParams.nFibers;
                    ymin = cParams.minyCoor;
                    ymax = cParams.maxyCoor;
                    k    = 2*pi*n/(ymax-ymin);
                    fH   = @(x) sin(k*(x2(x)-ymin)+pi/2);
                    
                    obj.fHandle = @(x) max(fV(x),fH(x));

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

                 case 'DiagonalBars'

                     % === Horitzontals (0°) ===
                    h = cParams.barWidth;
                    yB2 = cParams.bottomBar_yMax;
                    yT1 = cParams.topBar_yMin;
                    yB1 = yB2 - h;
                    yT2 = yT1 + h;
                
                    fH1 = @(x) max(yB1 - x2(x), x2(x) - yB2);
                    fH2 = @(x) max(yT1 - x2(x), x2(x) - yT2);
                
                    % === Diagonals (+45° respecte a X) ===
                    % Factor per mantenir el gruix correcte
                    scale = sqrt(2);
                
                    % Distància vertical entre les horitzontals
                    deltaY = yT1 - yB2;
                
                    % Separació equivalent entre diagonals (en coordenades c)
                    % La relació entre la distància vertical i el desplaçament c és √2
                    deltaC = deltaY * sqrt(2);
                
                    % Tres diagonals simètriques: inferior, central, superior
                    c_center = 0.0;
                    c_bottom = c_center - deltaC;
                    c_top    = c_center + deltaC;
                
                    % Definició de les tres barres diagonals
                    fD1 = @(x) abs(x2(x) - x1(x) - c_bottom) - (h/2)*scale;
                    fD2 = @(x) abs(x2(x) - x1(x) - c_center) - (h/2)*scale;
                    fD3 = @(x) abs(x2(x) - x1(x) - c_top) - (h/2)*scale;
                
                    % === Combina totes ===
                    fH = @(x) min(min(fH1(x), fH2(x)), min(min(fD1(x), fD2(x)), fD3(x)));
                
                    obj.fHandle = fH;

                case 'DiagonalBarsv2'
                    % === Horitzontals (0°) ===
                    L=1;
                    h = cParams.barWidth;
                    yB2 = cParams.bottomBar_yMax;
                    yT1 = cParams.topBar_yMin;
                    yB1 = yB2 - h;
                    yT2 = yT1 + h;
                
                    fH1 = @(x) max(yB1 - x2(x), x2(x) - yB2);
                    fH2 = @(x) max(yT1 - x2(x), x2(x) - yT2);

                    % Distància vertical entre les horitzontals
                    deltaY = yT1 - yB2;

                    % === Diagonals (+45° respecte a X) ===
                    scale = sqrt(2); % Factor per mantenir el gruix correcte
                    deltaC = deltaY * sqrt(2); % La relació entre la distància vertical i el desplaçament c és √2

                    % Tres diagonals simètriques: inferior, central, superior
                    c_center = 0.0;
                    c_bottom = c_center - L/2;
                    c_top    = c_center + L/2;
                
                    % Definició de les tres barres diagonals
                    fD1 = @(x) abs(x2(x) - x1(x) - c_bottom) - (h/2)*scale;
                    fD2 = @(x) abs(x2(x) - x1(x) - c_center) - (h/2)*scale;
                    fD3 = @(x) abs(x2(x) - x1(x) - c_top) - (h/2)*scale;


                    fH = @(x) min(min(fH1(x), fH2(x)), min(min(fD1(x), fD2(x)), fD3(x)));

                    obj.fHandle = fH;
                   

               
                case 'CrossDiagonalBars' % -45°/+45°
                    h     = cParams.barWidth;
                    scale = sqrt(2);
                    L=1;
                
                    % +45°
                    c_center_p  = 0.0;
                    c_bottom_p  = c_center_p - L/2;
                    c_top_p     = c_center_p + L/2;
                
                    % -45°
                    c_center_m   = L;
                    c_bottom_m   = c_center_m - L/2;
                    c_top_m      = c_center_m + L/2;
                

                    % +45° diagonals
                    fP1 = @(x) abs(x2(x) - x1(x) - c_bottom_p) - (h/2)*scale;
                    fP2 = @(x) abs(x2(x) - x1(x) - c_center_p) - (h/2)*scale;
                    fP3 = @(x) abs(x2(x) - x1(x) - c_top_p)    - (h/2)*scale;
                
                    % -45° diagonals
                    fN1 = @(x) abs(x2(x) + x1(x) - c_bottom_m) - (h/2)*scale;
                    fN2 = @(x) abs(x2(x) + x1(x) - c_center_m) - (h/2)*scale;
                    fN3 = @(x) abs(x2(x) + x1(x) - c_top_m)    - (h/2)*scale;


                    fH = @(x) min(min(min(fP1(x), fP2(x)), fP3(x)), min(min(fN1(x), fN2(x)), fN3(x)) );
                    obj.fHandle = fH;

                    
                
                case 'DiagonalNFibers'
                    n    = cParams.nFibers;
                    xmin = cParams.minxCoor;
                    xmax = cParams.maxxCoor;
                    ymin = cParams.minyCoor;
                    ymax = cParams.maxyCoor;
                
                    % Longitud efectiva en direcció perpendicular a les diagonals
                    % (combinació de les dues coordenades)
                    L = (xmax - xmin + ymax - ymin) / sqrt(2);
                
                    % Definim el nombre d’ones equivalent al nombre de fibres
                    k = 2*pi*n / L;
                
                    % Definició del nivell de camp per línies a 45° (y = x + c)
                    % Aquí x2 - x1 és constant per a les línies de pendent 1
                    fH = @(x) sin(k * ((x2(x) - x1(x)) / sqrt(2)) + pi/2);
                
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

                case 'HorizontalNFibersV2'
                    h=0.1;
                    n    = cParams.nFibers;
                    ymin = cParams.minyCoor;
                    ymax = cParams.maxyCoor;
                    
                    % Distancia entre les rectes
                    dy= ymax/(n+1);

                    % Centres de cada recta
                    yCenters = dy*(1:n);
                    
                    % Crea el handle de totes les rectes
                    fList = cell(1,n);
                     for i = 1:n
                        y0 = yCenters(i);
                        fList{i} = @(x) abs(x2(x) - y0) - h/2;
                     end
                    
                     % Combinació de totes elles:
                    fH = @(x) fList{1}(x);
                    for i = 2:n
                        f_prev = fH;
                        fH = @(x) min(f_prev(x), fList{i}(x));
                    end
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

                case 'Empty'
                    fH = @(x) 1*ones(size(x1(x)));
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
