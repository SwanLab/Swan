classdef SuperEllipseExponentPlotter < handle
    
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
        
        function obj = SuperEllipseExponentPlotter(cParams)
            obj.init(cParams);            
        end
        
        function plot(obj)
            obj.plotXiRhoPlots();
            obj.plotMxMyPlots();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = cParams.fileName;
            obj.title    = cParams.title;
            obj.rhoV = cParams.rhoV;
            obj.xiV  = cParams.xiV;
            obj.qMean = cParams.qMean;            
            obj.value = cParams.value;            
        end
        
        function plotXiRhoPlots(obj)
            s.fileName = obj.fileName;
            s.title = obj.title;                        
            s.rhoV = obj.rhoV;
            s.xiV  = obj.xiV;
            s.value = obj.value;            
            p =  SuperEllipseExponentTxiRhoPlotter(s);
            p.plot();
        end
        
        function plotMxMyPlots(obj)          
            s.fileName = obj.fileName;
            s.title = obj.title;                        
            s.rhoV = obj.rhoV;
            s.xiV  = obj.xiV;
            s.qMean = obj.qMean;            
            s.value = obj.value;                        
            p =  SuperEllipseExponentMxMyPlotter(s);
            p.plot();            
        end
        
    end
    
end