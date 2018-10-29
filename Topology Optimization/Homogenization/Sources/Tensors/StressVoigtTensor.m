classdef StressVoigtTensor < VoigtTensor
    
    properties (Access = private)
     value
    end
    
    methods (Access = public)
        
        function obj = StressVoigtTensor()
            
        end
        
        function v = getValue(obj)
            v = obj.value;           
        end
        
        function setValue(obj,v)
            obj.value = v;            
        end
        
    end
    
end

