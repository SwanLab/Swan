classdef RHSintegratorFactory < handle

   methods (Access = public, Static)
       
       function obj = create(cParams)
           switch cParams.type
               case 'ShapeFunction'
                   % Computes the RHS using the NODAL FORCES and SHAPE
                   % FUNCTIONS
                   obj = RHSintegrator_ShapeFunction(cParams);
               case 'CutMesh'
                   % Computes the RHS using the NODAL FORCES and SHAPE
                   % FUNCTIONS for CUT meshes
                   obj = RHSintegrator_CutMesh(cParams);
               case 'ShapeDerivative'
                   % Computes the RHS using the NODAL FORCES and SHAPE
                   % FUNCTIONS' DERIVATIVE
                   obj = RHSintegrator_ShapeDerivative(cParams);
               case 'Elastic'
                   % Computes the RHS for ELASTIC problems
                   switch cParams.scale
                       case 'MACRO'
                            obj = RHSintegrator_ElasticMacro(cParams);
                       case 'MICRO'
                            obj = RHSintegrator_ElasticMicro(cParams);
                   end
           end
       end
       
   end
    
end