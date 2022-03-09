classdef LHSintegratorFactory < handle

   methods (Access = public, Static)
       
       function obj = create(cParams)
           switch cParams.type
               case 'MassMatrix'
                   %(computeElementalLHS) N*N
                   %+ assembleMatrix (LHSintegrator through integrator)
                   obj = LHSintergrator_Mass(cParams);
               case 'StiffnessMatrix'
                   %(computeElementalLHS) dN*dN
                   %+ assembleMatrix (LHSintegrator through integrator)
                   obj = LHSintergrator_Stiffness(cParams);
               case 'ElasticStiffnessMatrix'
                   %(computeElementalLHS) dN*C*dN (B*C*B)
                   %+ assembleMatrix (LHSintegrator through integrator)
                   obj = LHSintergrator_StiffnessElastic(cParams);
               case 'ElasticStiffnessMatrixOld'
                   % elemntal B + assamly --> globalB
                   % elemntal C + assamly --> globalC
                   % global B'*C*B
                   obj = LHSintergrator_StiffnessElasticStoredB(cParams);
                   %globalB in contructor
               case 'AnisotropicStiffnessMatrix'
                   % dB'*Celas*dB
                   obj = LHSintergratorAnisotropicStiffness(cParams);
           end
           
       end
       
   end
    
end