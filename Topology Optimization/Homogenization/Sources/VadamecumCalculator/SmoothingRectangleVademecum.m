classdef SmoothingRectangleVademecum < handle
    
    properties (Access = private)
        outPutPath
        smoothVad
        nonSmoothVad
        difDB
    end
    
    methods (Access = public)
        
        function obj = SmoothingRectangleVademecum()
            obj.init();           
            obj.smoothVad    = obj.computeVademecum('SmoothRectangle');
            obj.nonSmoothVad = obj.computeVademecum('Rectangle');
            obj.postprocessDifVademecum();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.createOutPutPath();
        end
        
        function createOutPutPath(obj)
            firstPart  = fullfile( '/home','alex','Dropbox');
            secondPart = fullfile('Amplificators','Images','MicroWithHole/');
            obj.outPutPath = fullfile(firstPart,secondPart);
        end        
        
        function vc = computeVademecum(obj,fileName)
            d.fileName = fileName;
            d.outPutPath = fullfile(obj.outPutPath,[d.fileName,'/']);
            vc = VademecumComputerAndPlotter(d);
            vc.compute();
        end
                
        function postprocessDifVademecum(obj)
            d.fileName   = 'ReactangleDifference';
            d.outPutPath = fullfile(obj.outPutPath,'Difference','/');
            d.smoothVD   = obj.smoothVad.vademecumData;
            d.smoothVD.feasibleIndex = obj.smoothVad.getFeasibleIndex();            
            d.nonSmoothVD = obj.nonSmoothVad.vademecumData;
            d.nonSmoothVD.feasibleIndex = obj.nonSmoothVad.getFeasibleIndex();                                    
            vc = VademecumDifferenceComputerAndPlotter(d);
            vc.compute();
        end             

    end
    
end