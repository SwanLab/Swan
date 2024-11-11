classdef DualVariable < handle
    
    properties (Access = public)
       fun
    end
    
    properties (Access = private)
       nConstraints
       valueOld
    end
    
   methods (Access = public, Static)
       
       function obj = DualVariable(cParams)
           obj.nConstraints = cParams.nConstraints;
           obj.fun.fValues  = zeros(obj.nConstraints,1);
       end
       
   end
   
   methods (Access = public)
      
       function restart(obj)
           obj.fun.fValues = obj.valueOld;
       end

       function update(obj,lVal)
           obj.fun.fValues = lVal;
       end
       
       function updateOld(obj)
           obj.valueOld = obj.fun.fValues;
       end
       
   end

end