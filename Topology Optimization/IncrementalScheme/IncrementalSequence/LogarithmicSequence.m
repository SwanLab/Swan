classdef LogarithmicSequence < IncrementalSequence
    
    methods (Access = protected)
        
        function generateAlphaSequence(obj)
             obj.alpha = logspace(obj.x0,obj.x1,obj.nSteps);
        end
        
    end
    
end