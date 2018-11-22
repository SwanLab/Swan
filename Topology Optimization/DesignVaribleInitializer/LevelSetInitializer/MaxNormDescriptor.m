classdef MaxNormDescriptor < NormDescriptor
      
    methods (Access = protected)
        
        function d = computeDistance(obj)
            x = obj.pos;
            d = max(abs(x),[],2) + 1e-14 - 1;
            obj.dist = d;
        end
        
    end
    
end

