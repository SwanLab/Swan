classdef  RectangleVademecumComputerForGivenVolume < ...
        VademecumComputerForGivenVolume
    
    methods (Access = public)
        
        function obj = RectangleVademecumComputerForGivenVolume(d)
            obj.init(d);
            obj.fileName = 'Rectangle';            
        end
        
    end
    
    methods (Access = protected)
        
        function findInclusionLengthForCertainVolume(obj)
            obj.my = obj.computeInclusionLengthForRectangle();
        end
        
    end
     
    
end