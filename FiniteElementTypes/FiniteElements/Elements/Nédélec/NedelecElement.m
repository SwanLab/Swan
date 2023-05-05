classdef NedelecElement < handle
    
    methods (Access = public, Static)
        
        function element = create(d)
            switch d
                case 1
                case 2
                    element = NedelecElement2D();
                case 3
                    element = NedelecElement3D();
            end
        end
    end
    
    methods (Access = public, Abstract)
        
    end
    
end