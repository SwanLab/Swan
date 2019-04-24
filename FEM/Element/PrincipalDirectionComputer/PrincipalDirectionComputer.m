classdef PrincipalDirectionComputer < handle
   
   properties (Access = public)
       direction 
       directionFunction       
   end
    
    
   methods (Access = public, Static)
       
       function obj = create(cParams)
           f = PrincipalDirectionComputerFactory();
           obj = f.create(cParams);           
       end       
       
   end
   
   methods (Access = public, Abstract)
      compute(obj) 
   end
   
    
end