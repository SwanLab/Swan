classdef DualVariable < handle
    
    properties (Access = public)
       value 
    end
    
    properties (Access = private)
       nConstraints 
       valueOld       
    end
    
   methods (Access = public, Static)
       
       function obj = DualVariable(cParams)
           obj.nConstraints = cParams.nConstraints;
           obj.value        = zeros(1,obj.nConstraints);           
       end
       
   end
   
   methods (Access = public)
      
       function restart(obj)
           obj.value = obj.valueOld;
       end
       
       function updateOld(obj)
           obj.valueOld = obj.value;
       end
       
   end

    
end