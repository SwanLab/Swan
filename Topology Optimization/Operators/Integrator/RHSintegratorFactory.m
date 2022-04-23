classdef RHSintegratorFactory < handle

   methods (Access = public, Static)
       
       function obj = create(cParams)
           switch cParams.type
               case 'ShapeFunction'
                   % Computes the RHS using the NODAL FORCES and SHAPE
                   % FUNCTIONS
                   obj = RHSintegrator_ShapeFunction(cParams);
               case 'ShapeFunctionCell'
                   % Computes the RHS using the NODAL FORCES and SHAPE
                   % FUNCTIONS for CUT meshes
                   obj = RHSintegrator_ShapeFunctionCutMesh(cParams);
               case 'ShapeDerivative'
                   % Computes the RHS using the NODAL FORCES and SHAPE
                   % FUNCTIONS' DERIVATIVE
                   obj = RHSintegrator_ShapeDerivative(cParams);
               case 'Elastic'
                   % Computes the RHS for ELASTIC problems
                   obj = RHSintegrator_Elastic(cParams);
           end
       end
       
   end
    
end