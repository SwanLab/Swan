classdef SuperEllipseExponentMxMyPlotter < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        fileName
        title
        rhoV
        xiV
        qMean
        value
    end
    
    methods (Access = public)
        
        function obj = SuperEllipseExponentMxMyPlotter(cParams)
            obj.init(cParams)
            
        end
        
         function plot(obj)
            obj.plotContour();                        
            obj.plotTriSurf();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = [cParams.fileName,'MxMy'];
            obj.title = cParams.title;                                    
            obj.rhoV = cParams.rhoV;
            obj.xiV  = cParams.xiV;
            obj.qMean = cParams.qMean;
            obj.value = cParams.value;               
        end
        
        function plotTriSurf(obj)
            s.fileName = obj.fileName;
            s.title = obj.title;                        
            s.axisAdder = MxMyAxisAdder();         
            s.view = [-16 53];            
            p =  SuperEllipseExponentTriSurfPlotter(s);
            x = SuperEllipseParamsRelator.mx(obj.xiV,obj.rhoV,obj.qMean);            
            y = SuperEllipseParamsRelator.my(obj.xiV,obj.rhoV,obj.qMean);      
            z = obj.value;                        
            p.plot(x,y,z);             
        end
        
        function plotContour(obj)
            s.fileName = obj.fileName;
            s.title = obj.title;            
            s.axisAdder = MxMyAxisAdder();
            p =  SuperEllipseExponentContourPlotter(s);
            x = SuperEllipseParamsRelator.mx(obj.xiV,obj.rhoV,obj.qMean);            
            y = SuperEllipseParamsRelator.my(obj.xiV,obj.rhoV,obj.qMean);      
            z = obj.value;                        
            p.plot(x,y,z);            
        end
        
    end
    
end