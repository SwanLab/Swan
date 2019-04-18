classdef InactiveConstraintsModifierNull < InactiveConstraintsModifier
    
    methods (Access = public)
        
        function modify(obj,cons,lambda,penalty)
           obj.constraint = cons; 
        end
        
    end
   
    
end
