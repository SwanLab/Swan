classdef MicroParams < DesignVariable
    
    methods (Access = public)
        
        function obj = MicroParams(cParams)
            obj.init(cParams);
        end
        
        function update(obj,value)
            obj.value = value;
        end
        
    end
    
    methods (Access = private)
        
    end
    
end
    
   