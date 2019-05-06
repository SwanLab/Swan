classdef MicroParams < DesignVariable
    
    methods (Access = public)
        
        function obj = MicroParams(cParams)
            obj.nVariables = 2;                        
            obj.init(cParams);
            obj.createValue();
        end
        
        function update(obj,value)
            obj.value = value;
        end
        
    end
    
    methods (Access = private)
        
        function createValue(obj)
            ndof = length(obj.mesh.coord(:,1));
            obj.value(:,1) = 0.02*ones(obj.nVariables*ndof,1);
        end
        
    end
    
end
    
   