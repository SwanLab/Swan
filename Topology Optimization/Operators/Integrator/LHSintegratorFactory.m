classdef LHSintegratorFactory < handle

   methods (Access = public, Static)
       
       function obj = create(cParams)
           switch cParams.type
               case 'MassMatrix'
                   obj = LHSintergrator_Mass(cParams);
               case 'StiffnessMatrix'
                   obj = LHSintergrator_Stiffness(cParams);
           end
           
       end       
       
   end
    
end