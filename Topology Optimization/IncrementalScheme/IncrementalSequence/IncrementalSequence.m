classdef IncrementalSequence < handle
    
    properties (GetAccess = public, SetAccess = protected)
        value
    end
    
    properties (Access = public)
        alpha
        initialValue
        finalValue
    end
    
    properties (Access = protected)
        x0
        x1
        nSteps
    end
    
    methods (Access = protected, Abstract)
        
        generateAlphaSequence(obj)
        
    end
    
    methods (Access = public)
        
        function update(obj,i)
            obj.value = (1-obj.alpha(i))*obj.initialValue + obj.alpha(i)*obj.finalValue;
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,x0,x1,nSteps,a0,a1)
            obj.x0 = x0;
            obj.x1 = x1;
            obj.nSteps = nSteps;
            obj.initialValue = a0;
            obj.finalValue = a1;
        end
        
    end
    
end

