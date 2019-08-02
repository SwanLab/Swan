classdef DesignVariable < handle
    
   properties (Access = public)
      value 
   end
   
   properties (GetAccess = public, SetAccess = private)
      xLength
   end
    
   properties (Access = private)
       
   end
    
   methods (Access = public)
       
       function obj = DesignVariable(cParams)
          obj.init(cParams)
       end    
       
   end
   
   methods (Access = private)
      
       function init(obj,cParams)
          obj.xLength = cParams.xLength;
          obj.value = zeros(obj.xLength,1); 
       end
       
   end
               
    
end