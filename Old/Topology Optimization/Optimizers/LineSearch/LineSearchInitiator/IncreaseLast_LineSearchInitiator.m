classdef IncreaseLast_LineSearchInitiator < LineSearchInitiator
    
    properties (Access = private)
        incrementFactor
        maxValue
        minValue
    end

    methods (Access = public)
        
        function obj = IncreaseLast_LineSearchInitiator(cParams)
            obj.init(cParams);
            obj.incrementFactor = cParams.incrementFactor;
            obj.maxValue        = cParams.maxValue;
            obj.minValue        = cParams.minValue;
        end
        
        function initStep = compute(obj,lastStep)
            if lastStep <= obj.minValue
                initStep = 5*lastStep;
            else
                k = obj.incrementFactor;
                initStep = k*lastStep;
                initStep = min(initStep,obj.maxValue);
            end
        end       
        
    end     
  
  
end

