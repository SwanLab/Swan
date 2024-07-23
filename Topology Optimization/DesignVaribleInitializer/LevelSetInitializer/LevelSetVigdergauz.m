classdef LevelSetVigdergauz < handle
    
    properties (Access = private)
        x
        y
        parameters
        levelSet
    end

    methods (Access = public)
        function obj = LevelSetVigdergauz(cParams)
            obj.rescaleCoordinates();
            obj.computeVigdergauzParameters(cParams);
            ls = obj.computeLevelSetInEllipticCoordinates();
            ls = obj.makeLevelSetNegativeOutRectangelEnvelope(ls);
            obj.levelSet = ls;
        end

        function fH = getFunctionHandle(obj)
            fH = obj.levelSet;
        end
    end
    
    methods (Access = private)

        function rescaleCoordinates(obj)
            obj.x = @(x) (x(1,:,:)-0.5);
            obj.y = @(x) (x(2,:,:)-0.5);
        end        
                
        function computeVigdergauzParameters(obj,cParams)
            s    = cParams.vigdergauzSettings;
            s.cx = 1;
            s.cy = 1;
            v    = VigdergauzParameters.create(s);      
            v.compute();
            obj.parameters = v.parameters;
        end
        
        function ls = computeLevelSetInEllipticCoordinates(obj)
            R       = obj.parameters.R;            
            [xe,ye] = obj.computeEllipticCoordinates();
            ls      = @(x) (1-xe(x).^2).*(1-ye(x).^2) - R;
        end       
        
        function ls = makeLevelSetNegativeOutRectangelEnvelope(obj,ls)
            out     = obj.isOutsideRectangleEnvelope();
            ls      = @(x) ls(x).*(~out(x))-abs(ls(x)).*out(x);              
        end
        
        function itIs = isOutsideRectangleEnvelope(obj)
            xv               = obj.x;
            yv               = obj.y;              
            isSmallerThanMx  = @(x) abs(xv(x)) <= obj.parameters.mx/2;
            isSmallerThanMy  = @(x) abs(yv(x)) <= obj.parameters.my/2;
            isInsideRectange = @(x) isSmallerThanMx(x) & isSmallerThanMy(x);  
            itIs             = @(x) ~isInsideRectange(x);            
        end
        
        function [xe,ye] = computeEllipticCoordinates(obj)
            xv    = obj.x;
            yv    = obj.y;            
            rx    = obj.parameters.rx;
            ry    = obj.parameters.ry;
            FxMax = obj.parameters.FxMax;
            FyMax = obj.parameters.FyMax;  
            mx    = obj.parameters.mx;
            my    = obj.parameters.my;  
            xe    = @(x) ellipj(xv(x)/(mx/2)*FxMax,rx);
            ye    = @(x) ellipj(yv(x)/(my/2)*FyMax,ry);            
        end
        
    end
    
end

