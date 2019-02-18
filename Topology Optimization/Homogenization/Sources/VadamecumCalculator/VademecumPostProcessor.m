classdef VademecumPostProcessor < handle
    
    properties (SetAccess = private, GetAccess = public)
       feasibleIndex               
    end
    
    properties (Access = private)        
       postData
       outPutPath
       fileName
       dBForPlotter
    end
    
    methods (Access = public)
        
        function obj = VademecumPostProcessor(d)
            obj.init(d);
        end
        
        function postprocess(obj)
            obj.makePlots();
            obj.printData();            
        end
               
    end
    
    methods (Access = private)
        
        function init(obj,d)
            obj.postData   = d.postData;  
            obj.outPutPath = d.outPutPath;
            obj.fileName   = d.fileName;
        end
        
        function makePlots(obj)
            obj.computeDataBaseForPlotter();
            obj.makeTxiRhoPlot();            
            obj.makeMxMyPlot();
        end
        
        function computeDataBaseForPlotter(obj)
            d.mxV        = obj.postData.mxV;
            d.myV        = obj.postData.myV;
            d.C          = obj.postData.Ctensor;
            d.invP       = obj.postData.PinvTensor;
            d.volume     = obj.postData.volume;
            d.hasToPrint = true;
            d.outPutPath = obj.outPutPath;
            d.microName  = obj.fileName;
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
        
        function printData(obj)
            d.outPutPath = [obj.outPutPath,obj.fileName];
            d.Ctensor  = obj.postData.Ctensor;
            d.Ptensor  = obj.postData.Ptensor;
            d.volume   = obj.postData.volume;
            p = VademecumDataPrinter(d);
            p.print();
        end        
        
    end
   
    
end   