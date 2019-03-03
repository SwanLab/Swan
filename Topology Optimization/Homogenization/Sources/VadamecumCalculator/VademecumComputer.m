classdef VademecumComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        vademecumData 
    end
    
    properties (Access = private)
        fileName
        outPutPath
        vademecumPlotter        
    end
    
    methods (Access = public)
        
        function obj = VademecumComputer(d)
            obj.init(d)            
        end
        
        function compute(obj)
            obj.calculateVademecum();
            obj.plotVademecum();
            obj.printData();                        
        end
        
        function fI = getFeasibleIndex(obj)
          vp = obj.vademecumPlotter;
          fI = vp.feasibleIndex;
        end
    end
    
    methods (Access = private)
        
        function init(obj,d)
            obj.fileName   = d.fileName;
            obj.outPutPath = d.outPutPath;            
        end
        
        function calculateVademecum(obj)
            d.fileName   = obj.fileName;
            d.outPutPath = fullfile(obj.outPutPath,[d.fileName,'/']);
            vc = VademecumCalculator(d);
            obj.vademecumData = vc.getData(); 
            %a  = load(obj.fileName);
            %obj.vademecumData = a.d;            
        end
        
        function plotVademecum(obj)
            d = obj.vademecumData;
            obj.vademecumPlotter = VademecumPlotters(d);
            obj.vademecumPlotter.plot();
        end
        
        function printData(obj)
            d.outPutPath = [obj.outPutPath,obj.fileName];
            d.Ctensor  = obj.vademecumData.postData.Ctensor;
            d.Ptensor  = obj.vademecumData.postData.Ptensor;
            d.volume   = obj.vademecumData.postData.volume;
            p = VademecumDataPrinter(d);
            p.print();
        end        
        
    end
    
end