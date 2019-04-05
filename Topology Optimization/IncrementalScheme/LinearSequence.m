classdef LinearSequence < IncrementalSequence
    
    methods (Access = protected)
        
        function generateAlphaSequence(obj)
             obj.alpha = linspace(obj.x0,obj.x1,obj.nSteps);
        end
        
    end
    
end