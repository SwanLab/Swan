classdef LevelSetWidthDescriptor < handle

    
    methods (Access = protected, Static)
        
        function w = computeWidth(m,x)
            wx = max(x) - min(x);
            w = 0.5*m*wx;
        end
        
    end
    
    
end

