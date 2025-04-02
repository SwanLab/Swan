classdef VademecumComputerForGivenVolumeFactory
        
    methods (Access = public, Static)
        
        function v = create(d)
            
            switch d.vademecumCase
                case 'SmoothRectangle'
                   v = SmoothRectangleVademecumComputerForGivenVolume(d);
                case 'Rectangle'
                   v = RectangleVademecumComputerForGivenVolume(d);
            end
            
        end
        
    end
            
end