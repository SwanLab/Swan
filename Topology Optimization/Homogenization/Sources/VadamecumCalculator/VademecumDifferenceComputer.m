classdef VademecumDifferenceComputer < handle
    
    properties (Access = private)
       smoothVD 
       nonSmoothVD
       difVD     
       fileName       
       outPutPath       
    end

    methods (Access = public)
        
        function obj = VademecumDifferenceComputer(d)
            obj.init(d)
        end
        
        function compute(obj)
            obj.calculateVademecum();
            obj.plotVademecum();
        end
    end
    
    methods (Access = private)
        
        function init(obj,d)
            obj.fileName    = d.fileName;
            obj.outPutPath  = d.outPutPath;  
            obj.smoothVD    = d.smoothVD;
            obj.nonSmoothVD = d.nonSmoothVD;            
        end        
        
        function calculateVademecum(obj)
            dS = obj.smoothVD.postData;
            dN = obj.nonSmoothVD.postData;
            dD.volume     = (dS.volume     - dN.volume);
            dD.Ctensor    = (dS.Ctensor    - dN.Ctensor);
            dD.Ptensor    = (dS.Ptensor    - dN.Ptensor);
            dD.PinvTensor = (dS.PinvTensor - dN.PinvTensor);
            dD.mxV        = dS.mxV;
            dD.myV        = dS.myV;
            obj.difVD = dD;          
        end
        
        function plotVademecum(obj)
            obj.makeMxMyPlot();
            obj.makeTxiRhoPlot();               
        end
        
        function makeTxiRhoPlot(obj)
            d.smoothDB    = obj.smoothVD.postData;
            d.iS          = obj.smoothVD.feasibleIndex;
            d.nonSmoothDB = obj.nonSmoothVD.postData;
            d.iN          = obj.nonSmoothVD.feasibleIndex;
            d.microName  = obj.fileName;
            d.outPutPath = obj.outPutPath;
            d.hasToPrint = true;
            vp = VademecumTxiRhoPlotterDiff(d);
            vp.plot();
        end
        
        
        function makeMxMyPlot(obj)
            d = obj.difVD;
            d.C          = obj.difVD.Ctensor;
            d.invP       = obj.difVD.PinvTensor;
            d.volume     = obj.difVD.volume;
            d.microName  = obj.fileName;
            d.outPutPath = obj.outPutPath;             
            d.hasToPrint = true;            
            p = VademecumMxMyPlotter(d);
            p.plot();               
        end
        
    end
    
    
end