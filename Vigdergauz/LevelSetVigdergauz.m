classdef LevelSetVigdergauz < LevelSetCreator
    
    properties (Access = private)
        x
        y
        ax
        ay
        cx
        cy
        parameters
        vigdergauzType
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
            obj.vigdergauzType     = cParams.vigdergauzType;
            obj.vigdergauzDataBase = cParams.vigdergauzDataBase;        
            obj.cx = 1; 
            obj.cy = 1;            
        end
        
        function rescaleCoordinates(obj)
            xv = obj.nodeCoord(:,1);
            yv = obj.nodeCoord(:,2);
            obj.x = (1 - (-1))*(xv-0.5);
            obj.y = (1 - (-1))*(yv-0.5);            
        end        
                
        function computeVigdergauzParameters(obj)
            s.volume = obj.vigdergauzDataBase.volumeMicro;
            s.phi    = atan(obj.vigdergauzDataBase.superEllipseRatio);
            s.cx     = obj.cx;
            s.cy     = obj.cy;
            v = VigdergauzParametersFromThetaAndPhi(s);      
            obj.parameters = v.parameters;
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

