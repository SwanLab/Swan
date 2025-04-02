classdef DesignImagVariable < handle
    
   properties (Access = public)
      value 
      valueOld
      valueOldOld
   end
   
   properties (GetAccess = public, SetAccess = private)
      xLength
   end
    
   properties (Access = private)
       
   end
    
   methods (Access = public)
       
       function obj = DesignImagVariable(cParams)
          obj.init(cParams)
       end    
       
       function update(obj)
           obj.valueOldOld = obj.valueOld;
           obj.valueOld = obj.value;
       end
       
   end
   
   methods (Access = private)
      
       function init(obj,cParams)
          obj.xLength = cParams.xLength;
          obj.valueOldOld = zeros(obj.xLength,1);                     
          obj.valueOld    = zeros(obj.xLength,1);           
          obj.value       = zeros(obj.xLength,1); 
       end
       
   end
    
end