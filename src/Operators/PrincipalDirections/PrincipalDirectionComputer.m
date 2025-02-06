classdef PrincipalDirectionComputer < handle
   
   properties (Access = public)
       direction 
       principalStress
   end
   
   properties (Access = protected)
       eigenComputer
       ndim
       eType
   end
    
   methods (Access = public, Static)
       
       function obj = create(cParams)
           f = PrincipalDirectionComputerFactory();
           obj = f.create(cParams);
       end
       
   end
   
   methods (Access = protected)
       
       function init(obj,cParams)
          obj.eType = cParams.eigenValueComputer.type;
          obj.createEigenValueAndVectorFunction();
       end
       
   end
   
   methods (Access = private)
       
       function createEigenValueAndVectorFunction(obj)
           s.ndim = obj.ndim;
           s.type = obj.eType;
           eV = EigenValueAndVectorComputer.create(s);
           obj.eigenComputer = eV;
       end
       
   end
   
   methods (Access = public, Abstract)
      compute(obj)
   end
    
end