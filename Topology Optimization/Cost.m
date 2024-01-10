classdef Cost < CC
    
    properties (Access = private)
        weights
    end
    
    properties (GetAccess = public, SetAccess = private)
        value0
    end
    
    methods (Access = public)
        
        function obj = Cost(cParams)
            obj.weights = cParams.weights;
            obj.init(cParams);
        end
        
        function c = computeNonNormalizedValue(obj)
            c = 0;
            for iSF = 1:length(obj.shapeFunctions)
                s0 = obj.shapeFunctions{iSF}.value0;
                s  = obj.shapeFunctions{iSF}.value;
                sR = s*s0;
                c = c + obj.weights(iSF)*sR;
            end
        end
        
    end
    
    methods (Access = protected)
        
        function updateFields(obj,iSF)
            newValue = obj.weights(iSF)*obj.shapeFunctions{iSF}.value;
            newGrad  = obj.weights(iSF)*obj.shapeFunctions{iSF}.gradient.fValues;
            if isempty (obj.value0)
                obj.value0 = newValue;
            end
            obj.value = obj.value + newValue/obj.value0;
            obj.gradient = obj.gradient + newGrad/obj.value0;
        end
   
    end
    
end
