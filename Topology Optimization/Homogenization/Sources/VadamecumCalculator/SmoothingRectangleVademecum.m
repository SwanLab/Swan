classdef SmoothingRectangleVademecum < handle
    
    properties (Access = private)
        outPutPath
    end
    
    
    methods (Access = public)
        
        function obj = SmoothingRectangleVademecum()
            obj.init()
            obj.computeNonSmoothRectangle();
            obj.computeSmoothRectangle();            
        end
        
    end
    
    
    methods (Access = private)
        
        function init(obj)
            obj.createOutPutPath();            
        end
        
        function computeSmoothRectangle(obj)
            d.fileName = 'SmoothRectangle';
            d.outPutPath = fullfile(obj.outPutPath,[d.fileName,'/']);            
            VademecumCalculator(d);
        end
        
        function computeNonSmoothRectangle(obj)
            d.fileName = 'Rectangle';
            d.outPutPath = fullfile(obj.outPutPath,[d.fileName,'/']);
            VademecumCalculator(d);
        end
        
        function createOutPutPath(obj)
            firstPart  = fullfile( '/home','alex','Dropbox');
            secondPart = fullfile('Amplificators','Images','MicroWithHole/');
            obj.outPutPath = fullfile(firstPart,secondPart);            
        end
        
    end
    
end