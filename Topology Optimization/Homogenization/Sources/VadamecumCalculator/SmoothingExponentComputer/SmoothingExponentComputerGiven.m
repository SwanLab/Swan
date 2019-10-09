classdef SmoothingExponentComputerGiven < SmoothingExponentComputer
    
   properties (Access = private)
       q
   end
    
   methods (Access = public)
       
       function obj = SmoothingExponentComputerGiven(cParams)
           obj.init(cParams)
       end       
       
   end

   methods (Access = protected)

       function computeExponent(obj) 
            obj.value = obj.q;
       end
       
   end   
   
   methods (Access = private)
       
       function init(obj,cParams)
          obj.q = cParams.q;
       end
   end   
 
end