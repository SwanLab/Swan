classdef CustomSequence < IncrementalSequence
    
    properties (Access = public)
        factor
    end
    
    methods (Access = public)
        
        function obj = CustomSequence(x0,x1,nSteps,initialValue,finalValue)
            obj@IncrementalSequence(x0,x1,nSteps,initialValue,finalValue);
        end
        
    end
    
    methods (Access = protected)
        
        function generateAlphaSequence(obj)
            if obj.nSteps < 2
                x = x1;
            else
                iSteps = 0:obj.nSteps-1;
                x = 1-(1-iSteps/(obj.nSteps-1)).^(obj.factor);
                x = (x1-x0)*x + x0;
            end
            obj.alpha = x;
        end
        
    end
    
end