classdef LevelSetVigdergauz < LevelSetCreator
    
    properties (Access = private)
        x
        y
        ax
        ay
        cx
        cy
        parameters
        vigdergauzDataBase        
    end
    
    methods (Access = public)
        
        function obj = LevelSetVigdergauz(cParams)
            obj.init(cParams);   
            obj.compute(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function computeLevelSet(obj)          
            obj.rescaleCoordinates();
            obj.computeVigdergauzParameters();                  
            ls = obj.computeLevelSetInEllipticCoordinates();
            ls = obj.makeLevelSetNegativeOutRectangelEnvelope(ls);
            obj.levelSet = ls;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.vigdergauzDataBase = cParams.vigdergauzDataBase;        
            obj.cx = 1; 
            obj.cy = 1;            
        end
        
        function rescaleCoordinates(obj)
            xv = obj.nodeCoord(:,1);
            yv = obj.nodeCoord(:,2);
%             obj.x = (1 - (-1))*(xv-0.5);
%             obj.y = (1 - (-1))*(yv-0.5);  
            obj.x = (xv-0.5);
            obj.y = (yv-0.5);
        end        
                
        function computeVigdergauzParameters(obj)
            s = obj.vigdergauzDataBase;
            s.cx     = obj.cx;
            s.cy     = obj.cy;
            v = VigdergauzParameters.create(s);      
            v.compute();
            obj.parameters = v.parameters;
        end
        
        function ls = computeLevelSetInEllipticCoordinates(obj)
            R  = obj.parameters.R;            
            [xe,ye] = obj.computeEllipticCoordinates();
            ls(:,1) = (1-xe.^2).*(1-ye.^2) - R;
        end       
        
        function ls = makeLevelSetNegativeOutRectangelEnvelope(obj,ls)
            out = obj.isOutsideRectangleEnvelope();
            ls(out) = -abs(ls(out));                
        end
        
        function itIs = isOutsideRectangleEnvelope(obj)
            xv = obj.x;
            yv = obj.y;              
            isSmallerThanMx  = abs(xv) <= obj.parameters.mx/2;
            isSmallerThanMy  = abs(yv) <= obj.parameters.my/2;
            isInsideRectange =  isSmallerThanMx & isSmallerThanMy;  
            itIs = ~isInsideRectange;            
        end
        
        function [xe,ye] = computeEllipticCoordinates(obj)
            xv = obj.x;
            yv = obj.y;            
            rx = obj.parameters.rx;
            ry = obj.parameters.ry;
            FxMax = obj.parameters.FxMax;
            FyMax = obj.parameters.FyMax;  
            mx = obj.parameters.mx;
            my = obj.parameters.my;  
            xe = ellipj(xv/(mx/2)*FxMax,rx);
            ye = ellipj(yv/(my/2)*FyMax,ry);            
        end
        
    end
    
end

