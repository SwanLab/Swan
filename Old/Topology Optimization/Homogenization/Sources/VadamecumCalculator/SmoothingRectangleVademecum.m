classdef SmoothingRectangleVademecum < handle
    
    properties (Access = private)
        outPutPath
        smoothVad
        nonSmoothVad
        difDB
        
        firstVademecum
        secondVademecum
    end
    
    methods (Access = public)
        
        function obj = SmoothingRectangleVademecum()
            obj.init();           
            obj.smoothVad    = obj.computeVademecum(obj.firstVademecum);
            obj.nonSmoothVad = obj.computeVademecum(obj.secondVademecum);
            obj.postprocessDifVademecum();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.createOutPutPath();
            obj.firstVademecum = 'SuperEllipseQ2';            
            %obj.firstVademecum = 'SuperEllipseQOptAnalytic';
            obj.secondVademecum = 'SuperEllipseQMax';            
        end
        
        function createOutPutPath(obj)
            firstPart  = fullfile('/home','alex');
            secondPart = fullfile('git-repos','MicroStructurePaper','Other');
            obj.outPutPath = fullfile(firstPart,secondPart);
        end        
        
        function vc = computeVademecum(obj,fileName)
            d.fileName = fileName;
            d.outPutPath = fullfile(obj.outPutPath,[d.fileName,'/']);
            vc = VademecumComputerAndPlotter(d);
            vc.compute();
        end
                
        function postprocessDifVademecum(obj)
            fV = obj.firstVademecum;
            sV = obj.secondVademecum;
            folderPath = fullfile(obj.outPutPath,['DifferenceFrom',fV,'to',sV],'/');
            mkdir(folderPath);
            d.fileName   = 'ReactangleDifference';
            d.outPutPath = folderPath;
            d.smoothVD   = obj.smoothVad.vademecumData;
            d.smoothVD.feasibleIndex = obj.smoothVad.getFeasibleIndex();            
            d.nonSmoothVD = obj.nonSmoothVad.vademecumData;
            d.nonSmoothVD.feasibleIndex = obj.nonSmoothVad.getFeasibleIndex();                                    
            vc = VademecumDifferenceComputerAndPlotter(d);
            vc.compute();
        end             

    end
    
end