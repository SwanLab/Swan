classdef EpsilonSequence < IncrementalSequence
    
    methods (Access = protected)
        
        function generateAlphaSequence(obj)
            frac = 2;
            kmax = ceil(log10(obj.x0/obj.x1)/log10(frac));
            obj.alpha = obj.initialValue./frac.^(1:kmax);
        end
        
    end
    
end