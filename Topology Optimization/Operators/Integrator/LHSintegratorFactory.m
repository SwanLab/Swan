classdef LHSintegratorFactory < handle

   methods (Access = public, Static)

       function obj = create(cParams)
           switch cParams.type
               case 'MassMatrix'
                   obj = LHSintegrator_Mass(cParams);
               case 'MassMatrixVect'
                   obj = LHSintegrator_Mass_Vect(cParams);
               case 'BoundaryMassMatrix'
                   obj = LHSintegrator_MassBoundary(cParams);
               case 'StiffnessMatrix'
                   obj = LHSintegrator_Stiffness(cParams);
               case 'StiffnessMatrixVect'
                   obj = LHSintegrator_Stiffness_Vect(cParams);
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
               case 'MassBoundaryMass'
                   obj = LHSintegratorMassBoundaryMass(cParams);
               case 'AnisotropicStiffnessMassBoundaryMass'
                   cParams.stiffType = 'AnisotropicStiffnessMatrix';
                   obj = LHSintegratorStiffnessMassBoundaryMass(cParams);
               case 'AnisotropicStiffnessMass'
                   cParams.stiffType = 'AnisotropicStiffnessMatrix';
                   obj = LHSintegratorStiffnessMass(cParams);
               case 'Stokes'
                   obj = LHSintegrator_Stokes(cParams); % Anem a aquesta classe perquè nosaltres estem amb Stokes
               case 'Laplacian'
                   obj = LHSintegrator_Laplacian(cParams); % Venim de LHSintegrator_Stokes, perquè volem calcular el Laplacià
               case 'WeakDivergence'
                   obj = LHSintegrator_WeakDivergence(cParams);
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
