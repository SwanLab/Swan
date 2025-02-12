classdef Iterator < handle
    
   properties (GetAccess = public, SetAccess = private)
      value 
      maxIter
   end
    
   methods (Access = public)
       
       function obj = Iterator(cParams)
           obj.maxIter = cParams.maxIter;
           obj.value = 1;
       end
       
       function update(obj)
           obj.value = obj.value + 1;
       end
       
       function hasNotFinished = hasNotFinished(obj)
           hasNotFinished = obj.value < obj.maxIter;
       end
       
   end
    
end