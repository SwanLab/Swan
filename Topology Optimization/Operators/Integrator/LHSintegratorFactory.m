classdef LHSintegratorFactory < handle

   methods (Access = public, Static)

       function obj = create(cParams)
           switch cParams.type
               case 'MassMatrix'
                   obj = LHSintegrator_Mass(cParams);
               case 'BoundaryMassMatrix'
                   obj = LHSintegrator_MassBoundary(cParams);
               case 'StiffnessMatrix'
                   obj = LHSintegrator_Stiffness(cParams);
               case 'ElasticStiffnessMatrix'
                   obj = LHSintegrator_StiffnessElastic(cParams);
               case 'AnisotropicStiffnessMatrix'
                   obj = LHSintegrator_AnisotropicStiffness(cParams);
               case 'StiffnessMassBoundaryMass'
                   cParams.stiffType = 'StiffnessMatrix';
                   obj = LHSintegratorStiffnessMassBoundaryMass(cParams);
               case 'StiffnessMass'
                   cParams.stiffType = 'StiffnessMatrix';
                   obj = LHSintegratorStiffnessMass(cParams);
               case 'AnisotropicStiffnessMassBoundaryMass'
                   cParams.stiffType = 'AnisotropicStiffnessMatrix';
                   obj = LHSintegratorStiffnessMassBoundaryMass(cParams);
               case 'AnisotropicStiffnessMass'
                   cParams.stiffType = 'AnisotropicStiffnessMatrix';
                   obj = LHSintegratorStiffnessMass(cParams);
               case 'Stokes'
                   obj = LHSintegrator_Stokes(cParams);
               case 'Laplacian'
                   obj = LHSintegrator_Laplacian(cParams);
               case 'WeakDivergence'
                   obj = LHSintegrator_WeakDivergence(cParams);
               case 'AdvectionMatrix'
                   obj = LHSintegrator_Advection(cParams);
                   
               case 'StiffnessMatrixWithFunction'
                   obj = LHSintegratorFunctionStiffness(cParams);
               case 'MassMatrixWithFunction'
                   obj = LHSintegratorFunctionMass(cParams);
           end
       end
   end
end
