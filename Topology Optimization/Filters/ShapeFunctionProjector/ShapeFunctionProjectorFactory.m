classdef ShapeFunctionProjectorFactory < handle
    
   methods (Access = public, Static)
       
       function obj = create(cParams)
           switch cParams.type
               case 'TRIANGLE'
                   obj = ShapeFunctionProjector_ForTriangles(cParams);
               otherwise
                   obj = ShapeFunctionProjector_General(cParams);
           end
       end
       
   end
    
end