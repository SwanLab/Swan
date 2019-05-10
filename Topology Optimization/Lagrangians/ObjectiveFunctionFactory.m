classdef ObjectiveFunctionFactory < handle
    
    
   methods (Access = public, Static)
       
       function obj = create(cParams)
           switch cParams.type
               case 'AugmentedLagrangian'
                   obj = AugmentedLagrangian(cParams);
               case 'Lagrangian'
                   obj = Lagrangian(cParams);
           end
           
       end       
       
   end
    
    
    
end