classdef LHSintegratorFactory < handle

   methods (Access = public, Static)
       
       function obj = create(cParams)
           switch cParams.type
               case 'MassMatrix'
                   % Computes the MASS MATRIX by first computing the
                   % elemental LHS (N*N) and then assembling the result
                   obj = LHSintegrator_Mass(cParams);
               case 'BoundaryMassMatrix'
                   % Integrates the mass matrix over the boundary elements
                   % of the mesh
                   obj = LHSintegrator_MassBoundary(cParams);
               case 'StiffnessMatrix'
                   % Computes the STIFFNESS MATRIX by first computing the
                   % elemental LHS (dN*dN) and then assembling the result
                   obj = LHSintegrator_Stiffness(cParams);
               case 'ElasticStiffnessMatrix'
                   % Computes the ELASTIC STIFFNESS MATRIX by first
                   % computing the elemental LHS (dN*C*dN / B*C*B) and then
                   % assembling the result
                   obj = LHSintegrator_StiffnessElastic(cParams);
               case 'ElasticStiffnessMatrixOld'
                   % Computes the ELASTIC STIFFNESS MATRIX by first
                   % computing the *global* B matrix, the *global* C
                   % matrix, and then multiplying B*C*B globally.
                   % Per results, it is less efficient.
                   obj = LHSintegrator_StiffnessElasticStoredB(cParams);
                   %globalB in contructor
               case 'AdvectionMatrix'
                   %cross(b,grad(b))
                   obj = LHSintegrator_Advection(cParams);
               case 'DiffReactRobin'
                   % Creates a composite LHS for DiffReact problems with
                   % the ROBIN TERM. Includes a stiffness matrix, a mass
                   % matrix, and a boundary mass matrix.
                   obj = LHSintegrator_DiffReactRobin(cParams);
               case 'DiffReactNeumann'
                   % Creates a composite LHS for DiffReact problems with
                   % NO ROBIN TERM. Includes a stiffness matrix and a mass
                   % matrix, and NO boundary mass matrix.
                   obj = LHSintegrator_DiffReactNeumann(cParams);
               case 'Stokes'
                   obj = LHSintegrator_Stokes(cParams);

               case 'Laplacian'
                   obj = LHSintegrator_Laplacian(cParams);

               case 'StokesD'
                   obj = LHSintegrator_StokesD(cParams);
           end
       end
       
   end
    
end
