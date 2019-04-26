classdef MicroParams < DesignVariable
    
    methods (Access = public)
        
        function obj = MicroParams(cParams)
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
            obj.value(:,1) = 0.5*ndof;
            %obj.value(:,2) = 0.5*ndof;
        end
        
    end
    
end
    
   