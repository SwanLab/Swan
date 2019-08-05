classdef SplittingAlgorithmFactory < handle
    
   methods (Access = public, Static)
      
       function obj = create(cParams)
           switch cParams.type
               case 'ForwardBackward'
                   obj = ForwardBackward(cParams);
               case 'AcceleratedForwardBackward'
                   obj = AcceleratedForwardBackward(cParams);
           end

       end       
       
   end
    
    
end