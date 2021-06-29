classdef VademecumPlotters < handle
    
    properties (SetAccess = private, GetAccess = public)
       feasibleIndex               
    end
    
    properties (Access = private)        
       vademecumData
       outPutPath
       fileName
       dBForPlotter
    end
    
    methods (Access = public)
        
        function obj = VademecumPlotters(d)
            obj.init(d);
        end
        
        function plot(obj)
            obj.computeDataBaseForPlotter();
            obj.makeTxiRhoPlot();            
            obj.makeMxMyPlot();
        end
               
    end
    
    methods (Access = private)
        
        function init(obj,d)
            obj.vademecumData = d;  
        end
                
        function computeDataBaseForPlotter(obj)
            d.mxV        = obj.vademecumData.mxV;
            d.myV        = obj.vademecumData.myV;
            d.C          = obj.vademecumData.C;
            d.invP       = obj.vademecumData.invP;
            d.volume     = obj.vademecumData.volume;
            d.hasToPrint = false;
            d.outPutPath = obj.vademecumData.outPutPath;
            d.microName  = obj.vademecumData.fileName;
            obj.dBForPlotter = d;            
        end
        
        function makeMxMyPlot(obj)
            d = obj.dBForPlotter;
            p = VademecumMxMyPlotter(d);
            p.plot(); 
        end
        
        function makeTxiRhoPlot(obj)
            d = obj.dBForPlotter;            
            p = VademecumTxiRhoPlotter(d);
            p.plot();            
            obj.feasibleIndex = p.feasibleIndex;            
        end
        
    end
   
    
end   