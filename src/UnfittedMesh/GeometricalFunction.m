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
            % ls        = aFun.project('P1');

            sF.trial = LagrangianFunction.create(m,1,'P1');
            sF.mesh = m;
            filter = FilterLump(sF);
            ls = filter.compute(aFun,10);

        end

        function fxV = evaluate(obj,xV,m)
            s.fHandle = obj.fHandle;
            s.ndimf   = 1;
            s.mesh    = m;
            aFun      = AnalyticalFunction(s);
            fxV       = aFun.evaluate(xV);
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

                case 'Naca'
                    fH = @(x) obj.createNacaHole(x1(x),x2(x),cParams);
                    obj.fHandle = fH;

                case 'NacaInterior'
                    fH = @(x) obj.createNaca2(x1(x),x2(x),cParams);
                    obj.fHandle = fH;

                case 'NacaHole'
                    s      = cParams;
                    s.type = 'NacaInterior';
                    obj.computeInclusion(s);
                case 'LevelSetTest'
                    fH = @(x) obj.createLSTest(x1(x),x2(x),cParams);
                    obj.fHandle = fH;

                case 'LevelSet1'
                    fH = @(x) obj.createLS1(x1(x),x2(x),cParams);
                    obj.fHandle = fH;

                case 'LevelSet2'
                    fH = @(x) obj.createLS2(x1(x),x2(x),cParams);
                    obj.fHandle = fH;

                case 'LevelSet3'
                    fH = @(x) obj.createLS3(x1(x),x2(x),cParams);
                    obj.fHandle = fH;

                case 'LevelSet4'
                    fH = @(x) obj.createLS4(x1(x),x2(x),cParams);
                    obj.fHandle = fH;
                case 'LevelSet5'
                    fH = @(x) obj.createLS4(x1(x),x2(x),cParams);
                    obj.fHandle = fH;
                case 'LevelSet6'
                    fH = @(x) obj.createLS4(x1(x),x2(x),cParams);
                    obj.fHandle = fH;
                case 'LevelSet7'
                    fH = @(x) obj.createLS4(x1(x),x2(x),cParams);
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

        function fV = createNacaHole(x,y,s)
            c   = s.chord;
            p   = s.p;
            m   = s.m;
            t   = s.t;
            AoA = deg2rad(s.AoA);
        
            x0     = s.xLE;
            y0     = s.yLE/c;
            offsetX  = (x - x0)/c;
            offsetY  = y/c - y0;

            xNaca    = offsetX.*cos(AoA) - offsetY.*sin(AoA);
            yNaca    = offsetX.*sin(AoA) + offsetY.*cos(AoA);
        
            yc   = (xNaca>=0 & xNaca<=p).*(m./p^2.*(2*p*xNaca-xNaca.^2))+...
                    (xNaca>p & xNaca<=1).*(m./(1-p)^2.*((1-2*p)+2*p*xNaca-xNaca.^2));
            yt   = (xNaca>=0 & xNaca<=1).*(5*t*(0.2969*sqrt(xNaca)-0.1260*xNaca-0.3516*xNaca.^2+0.2843*xNaca.^3-0.1036*xNaca.^4));
            dydx = (xNaca>=0 & xNaca<=p).*(2*m/p^2.*(p-xNaca))+...
                    (xNaca>p & xNaca<=1).*(2*m/(1-p)^2.*(p-xNaca));
        
            theta = atan(dydx);
            yu    = yc + yt.*cos(theta);
            yl    = yc - yt.*cos(theta);
            
            f(:,:,:,1)   = yl - yNaca;
            f(:,:,:,2)   = yNaca - yu;
            f(:,:,:,3)   = xNaca - 1; 
            f(:,:,:,4)   = -xNaca;

            % No se quina d'aquestes dues funcionaria. Jo primer resoldria
            % el tema de la rotació i després revisem els detallets dels
            % forats petits

            %f  = f./max(abs(f));

            % f(:,:,:,1) = f(:,:,:,1)./max(abs(f(:,:,:,1)),[],'all');
            % f(:,:,:,2) = f(:,:,:,2)./max(abs(f(:,:,:,2)),[],'all');
            % f(:,:,:,3) = f(:,:,:,3)./max(abs(f(:,:,:,3)),[],'all');
            % f(:,:,:,4) = f(:,:,:,4)./max(abs(f(:,:,:,4)),[],'all');

            
            fV = -max(f,[],4);       
        end

        


        function fV = createNaca2(x,y,s)

            % STEP 1. Define Naca level-sets
            c   = s.chord;
            p   = s.p;
            m   = s.m;
            t   = s.t;
            AoA = deg2rad(s.AoA);
        
            x0     = s.xLE;
            y0     = s.yLE/c;
            offsetX  = (x - x0)/c;
            offsetY  = y/c - y0;

            xNaca    = offsetX.*cos(AoA) - offsetY.*sin(AoA);
            yNaca    = offsetX.*sin(AoA) + offsetY.*cos(AoA);
        
            yc   = (xNaca>=0 & xNaca<=p).*(m./p^2.*(2*p*xNaca-xNaca.^2))+...
                    (xNaca>p & xNaca<=1).*(m./(1-p)^2.*((1-2*p)+2*p*xNaca-xNaca.^2));
            yt   = (xNaca>=0 & xNaca<=1).*(5*t*(0.2969*sqrt(xNaca)-0.1260*xNaca-0.3516*xNaca.^2+0.2843*xNaca.^3-0.1036*xNaca.^4));
            dydx = (xNaca>=0 & xNaca<=p).*(2*m/p^2.*(p-xNaca))+...
                    (xNaca>p & xNaca<=1).*(2*m/(1-p)^2.*(p-xNaca));
        
            theta = atan(dydx);
            yu    = yc + yt.*cos(theta);
            yl    = yc - yt.*cos(theta);
            
            f(:,:,:,1)   = yl - yNaca;
            f(:,:,:,2)   = yNaca - yu;
            f(:,:,:,3)   = xNaca - 1; 
            f(:,:,:,4)   = -xNaca;
            
            % STEP 2. Make +/- signs more sparse
            f(:,:,:,1)   = -50.^(-3*f(:,:,:,1)) + 1;
            f(:,:,:,2)   = -50.^(-3*f(:,:,:,2)) + 1;
            f(:,:,:,3)   = -50.^(-18*f(:,:,:,3)) + 1;
            f(:,:,:,4)   = -50 .^(-18*f(:,:,:,4)) + 1;
            % f(:,:,:,1)   = -exp(-10*f(:,:,:,1)) + 1;
            % f(:,:,:,2)   = -exp(-10*f(:,:,:,2)) + 1;
            % f(:,:,:,3)   = -exp(-10*f(:,:,:,3)) + 1;
            % f(:,:,:,4)   = -exp(-10*f(:,:,:,4)) + 1;
            % 
            % STEP 3. Smooth the max function
            % ...
            % minF = min(f(:));
            % f = f - minF; 
            % 
            % p = 9;
            % 
            % fV = (f(:,:,:,1).^p + f(:,:,:,2).^p + f(:,:,:,3).^p + f(:,:,:,4).^p).^(1/p) + minF;
            

            % alpha = 200;
            % 
            % fVUp   = 0;
            % fVDown = 0;
            % 
            % for i = 1:4
            % 
            %     fun      = f(:,:,:,i);
            %     fVUp   = fVUp + fun.*exp(alpha*fun);
            %     fVDown = fVDown + exp(alpha*fun); 
            % 
            % end

            % fV = fVUp./fVDown;

            alpha = 2;
            fV = 1/alpha*log(exp(alpha*f(:,:,:,1)) + exp(alpha*f(:,:,:,2)) + exp(alpha*f(:,:,:,3)) + exp(alpha*f(:,:,:,4)));


            % alpha = 200;
            % fV = 1/alpha*log(1/4*(exp(alpha*f(:,:,:,1)) + exp(alpha*f(:,:,:,2)) + exp(alpha*f(:,:,:,3)) + exp(alpha*f(:,:,:,4))));


            % fV = max(f,[],4);     %Here no "-" sign is needed since we differentiate the inner naca from the hole naca

            % Crear el filtro gaussiano (Neteja tot darrere, la part separada)
            % sigma = 3;  % Ajustar según la necesidad
            % kernelSize = ceil(6 * sigma);  
            % if mod(kernelSize,2) == 0 
            %     kernelSize = kernelSize + 1;
            % end
            % 
            % 
            % kernel = fspecial('gaussian', [kernelSize kernelSize], sigma);
            % fV = imfilter(fV, kernel, 'same');

            % Detectar borde de salida
            % edgeMask = xNaca > 0.95 & xNaca < 1.05;  
            % 
            % % % Suavizar solo en el borde de salida
            % % fV(edgeMask) = conv2(fV(edgeMask), ones(1,5)/5, 'same');
            % 
            % % Levantar un poco más el borde de salida añadiendo un offset
            % offset = -0.05;  % Ajusta este valor según lo necesites
            % fV = fV + offset;
            % %fV(edgeMask) = fV(edgeMask) + offset;

        end

        function fV = createLSTest(x,y,s)

            c   = s.chord;
            p   = s.p;
            m   = s.m;
            t   = s.t;
            AoA = deg2rad(s.AoA);
        
            x0     = s.xLE;
            y0     = s.yLE/c;
            offsetX  = (x - x0)/c;
            offsetY  = y/c - y0;

            xNaca    = offsetX.*cos(AoA) - offsetY.*sin(AoA);
            yNaca    = offsetX.*sin(AoA) + offsetY.*cos(AoA);

            yc   = (xNaca>=0 & xNaca<=p).*(m./p^2.*(2*p*xNaca-xNaca.^2))+...
                    (xNaca>p & xNaca<=1).*(m./(1-p)^2.*((1-2*p)+2*p*xNaca-xNaca.^2));
            yt   = (xNaca>=0 & xNaca<=1).*(5*t*(0.2969*sqrt(xNaca)-0.1260*xNaca-0.3516*xNaca.^2+0.2843*xNaca.^3-0.1036*xNaca.^4));
            dydx = (xNaca>=0 & xNaca<=p).*(2*m/p^2.*(p-xNaca))+...
                    (xNaca>p & xNaca<=1).*(2*m/(1-p)^2.*(p-xNaca));


            theta = atan(dydx);
            yu    = yc + yt.*cos(theta);
            yl    = yc - yt.*cos(theta);
            
            %Primer intent fallat: amb producte sembla que no millora la
            %cosa, crec que això està fent lo mateix com maximització,
            %esforçant el cumpliment de punts
            fV   = (yNaca - yu).*(yl - yNaca);
            % 
            % f(:,:,:,1)   = yl - yNaca;
            % f(:,:,:,2)   = yNaca - yu;
            % f(:,:,:,3)   = xNaca - 1; 
            % f(:,:,:,4)   = -xNaca;

            % f(:,:,:,1)   = -exp(-50*f(:,:,:,1)) + 1;
            % f(:,:,:,2)   = -exp(-50*f(:,:,:,2)) + 1;
            % f(:,:,:,3)   = -exp(-50*f(:,:,:,3)) + 1;
            % f(:,:,:,4)   = -exp(-50*f(:,:,:,4)) + 1;

            % fV = -((yNaca - yc).^2 - yt.^2.*cos(theta).^2);

            %fV = -max(f,[],4);

        end



         function fV = createLS1(x,y,s)

            c   = s.chord;
            p   = s.p;
            m   = s.m;
            t   = s.t;
            AoA = deg2rad(s.AoA);
        
            x0     = s.xLE;
            y0     = s.yLE/c;
            offsetX  = (x - x0)/c;
            offsetY  = y/c - y0;

            xNaca    = offsetX.*cos(AoA) - offsetY.*sin(AoA);
            yNaca    = offsetX.*sin(AoA) + offsetY.*cos(AoA);
        
            yc   = (xNaca>=0 & xNaca<=p).*(m./p^2.*(2*p*xNaca-xNaca.^2))+...
                    (xNaca>p & xNaca<=1).*(m./(1-p)^2.*((1-2*p)+2*p*xNaca-xNaca.^2));
            yt   = (xNaca>=0 & xNaca<=1).*(5*t*(0.2969*sqrt(xNaca)-0.1260*xNaca-0.3516*xNaca.^2+0.2843*xNaca.^3-0.1036*xNaca.^4));
            dydx = (xNaca>=0 & xNaca<=p).*(2*m/p^2.*(p-xNaca))+...
                    (xNaca>p & xNaca<=1).*(2*m/(1-p)^2.*(p-xNaca));
        
            theta = atan(dydx);
            yl    = yc - yt.*cos(theta);
            
             fV   = -(yl - yNaca);

            % f(:,:,:,1)   = -(yl - yNaca);
            % f(:,:,:,2)   = -xNaca;
            % fV           = max(f,[],4);
            %fV           = -fV;

         end

             function fV = createLS2(x,y,s)

            c   = s.chord;
            p   = s.p;
            m   = s.m;
            t   = s.t;
            AoA = deg2rad(s.AoA);
        
            x0     = s.xLE;
            y0     = s.yLE/c;
            offsetX  = (x - x0)/c;
            offsetY  = y/c - y0;

            xNaca    = offsetX.*cos(AoA) - offsetY.*sin(AoA);
            yNaca    = offsetX.*sin(AoA) + offsetY.*cos(AoA);
        
            yc   = (xNaca>=0 & xNaca<=p).*(m./p^2.*(2*p*xNaca-xNaca.^2))+...
                    (xNaca>p & xNaca<=1).*(m./(1-p)^2.*((1-2*p)+2*p*xNaca-xNaca.^2));
            yt   = (xNaca>=0 & xNaca<=1).*(5*t*(0.2969*sqrt(xNaca)-0.1260*xNaca-0.3516*xNaca.^2+0.2843*xNaca.^3-0.1036*xNaca.^4));
            dydx = (xNaca>=0 & xNaca<=p).*(2*m/p^2.*(p-xNaca))+...
                    (xNaca>p & xNaca<=1).*(2*m/(1-p)^2.*(p-xNaca));
        
            theta = atan(dydx);
            yu    = yc + yt.*cos(theta);
            
            fV   = -(yNaca - yu);

            % f(:,:,:,1)   = yNaca - yu;
            % f(:,:,:,2)   = -xNaca;
            % fV           = max(f,[],4);
            % fV           = - fV;
        end


        function fV = createLS3(x,y,s)

            c   = s.chord;
            p   = s.p;
            m   = s.m;
            t   = s.t;
            AoA = deg2rad(s.AoA);
        
            x0     = s.xLE;
            y0     = s.yLE/c;
            offsetX  = (x - x0)/c;
            offsetY  = y/c - y0;

            xNaca    = offsetX.*cos(AoA) - offsetY.*sin(AoA);
            yNaca    = offsetX.*sin(AoA) + offsetY.*cos(AoA);
        
            yc   = (xNaca>=0 & xNaca<=p).*(m./p^2.*(2*p*xNaca-xNaca.^2))+...
                            (xNaca>p & xNaca<=1).*(m./(1-p)^2.*((1-2*p)+2*p*xNaca-xNaca.^2));
            yt   = (xNaca>=0 & xNaca<=1).*(5*t*(0.2969*sqrt(xNaca)-0.1260*xNaca-0.3516*xNaca.^2+0.2843*xNaca.^3-0.1036*xNaca.^4));
            dydx = (xNaca>=0 & xNaca<=p).*(2*m/p^2.*(p-xNaca))+...
                   (xNaca>p & xNaca<=1).*(2*m/(1-p)^2.*(p-xNaca));
        
            theta = atan(dydx);
            yu    = yc + yt.*cos(theta);
            fV   = -(yNaca - yu);



        end

        
        function fV = createLS4(x,y,s)

            c   = s.chord;
            AoA = deg2rad(s.AoA);
        
            x0     = s.xLE;
            y0     = s.yLE/c;
            offsetX  = (x - x0)/c;
            offsetY  = y/c - y0;

            xNaca    = offsetX.*cos(AoA) - offsetY.*sin(AoA);
            
            fV   = xNaca - 1; 

        end

        function fV = createLS5(x,y,s)

            c   = s.chord;
            AoA = deg2rad(s.AoA);
        
            x0     = s.xLE;
            y0     = s.yLE/c;
            offsetX  = (x - x0)/c;
            offsetY  = y/c - y0;

            xNaca    = offsetX.*cos(AoA) - offsetY.*sin(AoA);
            
            fV   = -xNaca;

        end

        function fV = createLS6(x,y,s)

            c   = s.chord;
            p   = s.p;
            m   = s.m;
  
            AoA = deg2rad(s.AoA);
        
            x0     = s.xLE;
            y0     = s.yLE/c;
            offsetX  = (x - x0)/c;
            offsetY  = y/c - y0;

            xNaca    = offsetX.*cos(AoA) - offsetY.*sin(AoA);
        
            yc   = (xNaca>=0 & xNaca<=p).*(m./p^2.*(2*p*xNaca-xNaca.^2))+...
                            (xNaca>p & xNaca<=1).*(m./(1-p)^2.*((1-2*p)+2*p*xNaca-xNaca.^2));
            fV   = yc;

        end

        function fV = createLS7(x,y,s)

            c   = s.chord;
            t   = s.t;
            AoA = deg2rad(s.AoA);
        
            x0     = s.xLE;
            y0     = s.yLE/c;
            offsetX  = (x - x0)/c;
            offsetY  = y/c - y0;

            xNaca    = offsetX.*cos(AoA) - offsetY.*sin(AoA);
            yt   = (xNaca>=0 & xNaca<=1).*(5*t*(0.2969*sqrt(xNaca)-0.1260*xNaca-0.3516*xNaca.^2+0.2843*xNaca.^3-0.1036*xNaca.^4));

        
            fV   = yt;

        end

        function fV = createLS8(x,y,s)

            c   = s.chord;
            p   = s.p;
            m   = s.m;
            t   = s.t;
            AoA = deg2rad(s.AoA);
        
            x0     = s.xLE;
            y0     = s.yLE/c;
            offsetX  = (x - x0)/c;
            offsetY  = y/c - y0;

            xNaca    = offsetX.*cos(AoA) - offsetY.*sin(AoA);
        
            yc   = (xNaca>=0 & xNaca<=p).*(m./p^2.*(2*p*xNaca-xNaca.^2))+...
                            (xNaca>p & xNaca<=1).*(m./(1-p)^2.*((1-2*p)+2*p*xNaca-xNaca.^2));
            yt   = (xNaca>=0 & xNaca<=1).*(5*t*(0.2969*sqrt(xNaca)-0.1260*xNaca-0.3516*xNaca.^2+0.2843*xNaca.^3-0.1036*xNaca.^4));
            dydx = (xNaca>=0 & xNaca<=p).*(2*m/p^2.*(p-xNaca))+...
                   (xNaca>p & xNaca<=1).*(2*m/(1-p)^2.*(p-xNaca));
        
            theta = atan(dydx);
            yu    = yc + yt.*cos(theta);
            fV   = yu;

        end
    
    end

end