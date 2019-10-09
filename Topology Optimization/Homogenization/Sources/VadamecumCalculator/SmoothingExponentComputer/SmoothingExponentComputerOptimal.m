classdef SmoothingExponentComputerOptimal < SmoothingExponentComputer
    
   properties (Access = private)
       alpha
       beta
       gamma   
       m1
       m2
   end
    
   methods (Access = public)
       
       function obj = SmoothingExponentComputerOptimal(cParams)
           obj.init(cParams)
       end       
       
   end

   methods (Access = protected)

       function computeExponent(obj) 
            a = obj.alpha;
            b = obj.beta;
            c = obj.gamma;           
            x = max(obj.m1,obj.m2);
            q = min(512,c*(1/(1-x^b))^a);  
            obj.value = q;
       end
       
   end   
   
   methods (Access = private)
       
       function init(obj,cParams)
          obj.m1 = cParams.m1;
          obj.m2 = cParams.m2;
          obj.alpha = 6;
          obj.beta  = 20;
          obj.gamma = 4;
       end
   end   
 
end