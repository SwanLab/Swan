classdef LinearSequence < IncrementalSequence
    
    methods (Access = public)
       
        function obj = LinearSequence(x0,x1,nSteps,initialValue,finalValue)
            obj.init(x0,x1,nSteps,initialValue,finalValue);
            obj.generateAlphaSequence();
        end                
        
    end
    
    methods (Access = protected)
        
        function generateAlphaSequence(obj)
             obj.alpha = linspace(obj.x0,obj.x1,obj.nSteps);
        end
        
    end
    
end