classdef LevelSetLipung < handle
    
    properties (Access = public)
        value
    end
    
    properties (Access = private)
        x
        y
        volume
        phi
        ax
        ay
        cx
        cy
        parameters
    end
    
    methods (Access = public)
        
        function obj = LevelSetLipung(cParams)
            obj.init(cParams)
            obj.compute();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.x      = cParams.x;
            obj.y      = cParams.y;
            obj.volume = cParams.volum;
            obj.phi    = cParams.phi;
            obj.cx = 1; 
            obj.cy = 1;            
        end
        
        function compute(obj)
            obj.rescaleCoordinates();
            obj.computeVigdergauzParameters();
            obj.computeLevelSet();
        end
        
        function computeVigdergauzParameters(obj)
            s.volume = obj.volume;
            s.phi    = obj.phi;
            s.cx     = obj.cx;
            s.cy     = obj.cy;
            v = VigdergauzParametersFromThetaAndPhi(s);      
            obj.parameters = v.parameters;
        end
        
        function rescaleCoordinates(obj)
            obj.x = (1 - (-1))*(obj.x-0.5);
            obj.y = (1 - (-1))*(obj.y-0.5);            
        end
         
        function computeLevelSet(obj)
            ls = obj.computeLevelSetInEllipticCoordinates();
            ls = obj.makeLevelSetNegativeOutRectangelEnvelope(ls);
            obj.value = ls;
        end
        
        function ls = computeLevelSetInEllipticCoordinates(obj)
            M  = obj.parameters.M;            
            [xe,ye] = obj.computeEllipticCoordinates();
            ls(:,1) = (1-xe.^2).*(1-ye.^2) - M;
        end       
        
        function ls = makeLevelSetNegativeOutRectangelEnvelope(obj,ls)
            out = obj.isOutsideRectangleEnvelope();
            ls(out) = -abs(ls(out));                
        end
        
        function itIs = isOutsideRectangleEnvelope(obj)
            xv = obj.x;
            yv = obj.y;              
            isSmallerThanRx  = abs(xv) <= obj.parameters.rx;
            isSmallerThanRy  = abs(yv) <= obj.parameters.ry;
            isInsideRectange =  isSmallerThanRx & isSmallerThanRy;  
            itIs = ~isInsideRectange;            
        end
        
        function [xe,ye] = computeEllipticCoordinates(obj)
            xv = obj.x;
            yv = obj.y;            
            mx = obj.parameters.mx;
            my = obj.parameters.my;
            FxMax = obj.parameters.FxMax;
            FyMax = obj.parameters.FxMax;  
            rx = obj.parameters.rx;
            ry = obj.parameters.ry;  
            xe = ellipj(xv/rx*FxMax,mx);
            ye = ellipj(yv/ry*FyMax,my);            
        end
        
    end
    
end

