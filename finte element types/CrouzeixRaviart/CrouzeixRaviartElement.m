classdef CrouzeixRaviartElement < handle
    
    methods (Access = public, Static)
        
        function element = create(d)
            switch d
                case 1
                    element = CrouzeixRaviart1D();
                case 2
                    element = CrouzeixRaviart2D();
            end
        end
    end
    
    methods (Access = public, Abstract)
        
    end
    
end