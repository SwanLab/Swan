classdef LHSIntegratorFactory < handle

   methods (Access = public, Static)

       function obj = create(cParams)
           switch cParams.type
               case 'MassMatrix'
                   obj = LHSIntegratorMass(cParams);
               case 'MassMatrixVect'
                   obj = LHSintegrator_Mass_Vect(cParams);
               case 'BoundaryMassMatrix'
                   obj = LHSIntegratorMassBoundary(cParams);
               case 'StiffnessMatrix'
                   obj = LHSIntegratorStiffness(cParams);
               case 'StiffnessMatrixVect'
                   obj = LHSintegrator_Stiffness_Vect(cParams);
               case 'ElasticStiffnessMatrix'
                   obj = LHSIntegratorStiffnessElastic(cParams);
               case 'AnisotropicStiffnessMatrix'
                   obj = LHSIntegratorAnisotropicStiffness(cParams);
               case 'StiffnessMassBoundaryMass'
                   cParams.stiffType = 'StiffnessMatrix';
                   obj = LHSIntegratorStiffnessMassBoundaryMass(cParams);
               case 'StiffnessMass'
                   cParams.stiffType = 'StiffnessMatrix';
                   obj = LHSIntegratorStiffnessMass(cParams);
               case 'MassBoundaryMass'
                   obj = LHSIntegratorMassBoundaryMass(cParams);
               case 'AnisotropicStiffnessMassBoundaryMass'
                   cParams.stiffType = 'AnisotropicStiffnessMatrix';
                   obj = LHSIntegratorStiffnessMassBoundaryMass(cParams);
               case 'AnisotropicStiffnessMass'
                   cParams.stiffType = 'AnisotropicStiffnessMatrix';
                   obj = LHSIntegratorStiffnessMass(cParams);
               case 'Stokes'
                   obj = LHSIntegratorStokes(cParams);
               case 'Laplacian'
                   obj = LHSIntegratorLaplacian(cParams);
               case 'WeakDivergence'
                   obj = LHSIntegratorWeakDivergence(cParams);
               case 'AdvectionMatrix'
                   obj = LHSintegrator_Advection(cParams);
                   
               case 'StiffnessMatrixWithFunction'
                   obj = LHSintegratorFunctionStiffness(cParams);
               case 'MassMatrixWithFunction'
                   obj = LHSintegratorFunctionMass(cParams);
               case 'MassMatrixWithFunctionDerivative'
                   obj = LHSintegratorFunctionDerivativeMass(cParams);                   
               case 'AdvectionMatrixWithFunctionDerivative'
                   obj = LHSintegratorFunctionDerivativeAdvection(cParams);  
               case 'AdvectionMatrixWithFunction'
                   obj = LHSintegratorFunctionAdvection(cParams); 
               case 'DivergenceMatrix'
                   obj = LHSintegratorDivergenceMatrix(cParams); 
           end
       end
   end
end
