classdef  RectangleVademecumComputerForGivenVolume < ...
        VademecumComputerForGivenVolume
    
    properties (Access = protected)
       qValue 
    end
    
    methods (Access = public)
        
        function obj = RectangleVademecumComputerForGivenVolume(d)
            obj.init(d);
            obj.fileName = 'Rectangle';            
            obj.qValue = Inf;
        end
        
    end
    
    methods (Access = protected)
        
        function findInclusionLengthForCertainVolume(obj)
            obj.my = obj.computeInclusionLengthForRectangle();
        end
        
    end
     
    
end