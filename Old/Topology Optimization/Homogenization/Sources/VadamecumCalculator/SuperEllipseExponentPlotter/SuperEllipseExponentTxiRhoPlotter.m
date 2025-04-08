classdef SuperEllipseExponentTxiRhoPlotter < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        fileName
        title
        rhoV
        xiV
        value   
    end
    
    methods (Access = public)
        
        function obj = SuperEllipseExponentTxiRhoPlotter(cParams)
            obj.init(cParams);            
        end
        
        function plot(obj)
            obj.plotContour();                        
            obj.plotTriSurf();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = [cParams.fileName,'XiRho'];
            obj.title = cParams.title;                                    
            obj.rhoV = cParams.rhoV;
            obj.xiV  = cParams.xiV;
            obj.value = cParams.value;            
        end
        
        function plotTriSurf(obj)
            s.fileName = obj.fileName;
            s.title = obj.title;                        
            s.axisAdder = XiRhoAxisAdder();   
            s.view = [-223 67];
            p =  SuperEllipseExponentTriSurfPlotter(s);
            x = obj.xiV;            
            y = obj.rhoV;
            z = obj.value;                        
            p.plot(x,y,z);                  
        end
        
        function plotContour(obj)
            s.fileName = obj.fileName;
            s.title = obj.title;            
            s.axisAdder = XiRhoAxisAdder(); 
            p =  SuperEllipseExponentContourPlotter(s);
            x = obj.xiV;            
            y = obj.rhoV;
            z = obj.value;                        
            p.plot(x,y,z);                  
        end
        
        function addAxis(obj)
            zlabel('q');            
        end
        
    end
    
end