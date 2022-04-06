classdef EntryIndex < handle
    
    properties (Access = public)
        values
    end
    
    properties (Access = private)
        NumberOfIndex
        InitialValue
        FinalValue
    end
    
    methods (Access = public)
        
        function  obj = EntryIndex(NumberOfIndex)
            obj.NumberOfIndex = NumberOfIndex;
            obj.init();
        end
        
        function init(obj)
            obj.InitialValue=1;
        end
        
        function compute(obj)
            obj.updateFinalValue()
            obj.computeValues()
            obj.updateInitialValue()
        end
    end
    
    methods (Access = private)
        
        function updateFinalValue(obj)
            obj.FinalValue = obj.InitialValue + obj.NumberOfIndex - 1;
        end
        
        function computeValues(obj)
            obj.values = obj.InitialValue:obj.FinalValue;
        end
        
        function updateInitialValue(obj)
            obj.InitialValue = obj.FinalValue + 1;
        end
        
    end
    
end

