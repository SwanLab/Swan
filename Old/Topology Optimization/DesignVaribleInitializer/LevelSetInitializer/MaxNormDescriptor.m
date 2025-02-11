classdef MaxNormDescriptor < NormDescriptor
    
    properties (Access = protected, Abstract)
       pos 
    end
      
    methods (Access = protected)
        
        function computeDistance(obj)
            x = obj.pos;
            d = max(abs(x),[],2);
            obj.dist = d;
        end
        
    end
    
end

