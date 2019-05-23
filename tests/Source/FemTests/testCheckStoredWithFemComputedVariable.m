classdef testCheckStoredWithFemComputedVariable <  testShowingError  ...
         & testFemComputation ...
         & testLoadStoredVariable ...
         & testStoredComputedChecker
         
     properties (Access = protected)
       tol = 1e-6;        
     end
end

