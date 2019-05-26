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
            %obj.value(:,1) = 0.87*ones(obj.nVariables*ndof,1);
            
           % m1 = 0.87*ones(ndof,1);
           % m2 = 0.5*ones(ndof,1);
            
            %V = 0.3;
            %r = 1;
            %f = sqrt((1-V)/r);
            %m1 = f*ones(ndof,1);
            %m2 = r*f*ones(ndof,1);                        
            
            m1 = 0.87*ones(ndof,1);
            m2 = 0.87*ones(ndof,1);                        
            
            
            obj.value(:,1) = [m1;m2];
            
            %obj.value(:,1) = 0.87*ones(obj.nVariables*ndof,1);
            %obj.value(:,1) = 0.02*ones(obj.nVariables*ndof,1);            
        end
        
    end
    
end
    
   