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
            ls        = aFun.project('toP1');
        end
    end

    methods (Access = private)
        function selectHandle(obj,cParams)
            switch cParams.type
                case 'Square'
                    l = cParams.length;
                    x0 = cParams.xCoorCenter;
                    y0 = cParams.yCoorCenter;
                    xm = x0-l/2;
                    xM = x0+l/2;
                    ym = y0-l/2;
                    yM = y0+l/2;
                    fH = @(x) (x(2,:,:)-ym).*(x(1,:,:)>=xm && x(1,:,:)<=xM) +...
                        (x(1,:,:)-xM).*(x(2,:,:)>=ym && x(2,:,:)<=yM) +...
                        (x(2,:,:)-yM).*(x(1,:,:)>=xm && x(1,:,:)<=xM) +...
                        (x(1,:,:)-xm).*(x(2,:,:)>=ym && x(2,:,:)<=yM);
                    obj.fHandle = fH;

                case 'Circle'
                    r  = cParams.radius;
                    x0 = cParams.xCoorCenter;
                    y0 = cParams.yCoorCenter;
                    fH = @(x) (x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2-r^2;
                    obj.fHandle = fH;
            end
        end
    end
end