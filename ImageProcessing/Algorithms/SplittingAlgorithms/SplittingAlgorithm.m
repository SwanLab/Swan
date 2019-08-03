classdef SplittingAlgorithm < handle
    
   methods (Access = public, Static)
      
       function obj = create(cParams)
          f = SplittingAlgorithmFactory();
          obj = f.create(cParams);
       end       
       
   end
    
    
end