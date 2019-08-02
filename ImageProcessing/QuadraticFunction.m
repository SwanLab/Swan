classdef QuadraticFunction < handle
    
   properties (GetAccess = public) 
      lipschitzConstant 
   end
    
   properties (Access = private)
      Av
      bv
      designVariable
      value0 
   end
    
   methods (Access = public)
       
      
       function obj = QuadraticFunction(cParams)
           obj.Av = cParams.A;
           obj.bv = cParams.b;
           obj.designVariable    = cParams.designVariable;           
           obj.computeAdimensionalValue();
           obj.computeLipschitzConstant(cParams);           
       end
       
       function c = computeCost(obj)
           x = obj.designVariable.value;
           A = obj.Av;
           b = obj.bv;
           r = A*x - b;
           c0 = obj.value0;
           c = 0.5*(r'*r)/c0;
       end
       
       function g = computeGradient(obj)
           x = obj.designVariable.value;
           A = obj.Av;
           b = obj.bv;
           r = A*x - b;
           c0 = obj.value0;
           g = (A'*r)/c0;
       end
       
       function computeAdimensionalValue(obj)
           b = obj.bv;
           obj.value0 = 0.5*(b'*b);
       end
       
       function computeLipschitzConstant(obj,cParams)
           L = cParams.lipschitzConstant;
           c0 = obj.value0;
           obj.lipschitzConstant = L/c0;
       end
       
   end
    
   
    
end