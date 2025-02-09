classdef FreeSequence < IncrementalSequence
    
    methods (Access = protected)
        
        function generateAlphaSequence(obj)
            x = zeros(1,obj.nSteps);
            x(end) = 1;
            obj.alpha = x;
        end
        
    end
    
end