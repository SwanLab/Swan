classdef testAcceleratedDenoisingEinstein <  testShowingError...
         & testImageProcessingComputation ...
         & testLoadStoredVariable ...
         & testStoredComputedChecker
    
    properties (Access = protected)
        testName = 'test_AcceleratedDenoisingEinstein';  
        variablesToStore = {'x'};
        tol = 5e-2;        
    end
       
    methods (Access = protected)
        
        function selectComputedVar(obj)
            obj.computedVar{1} = obj.imageProblem.optimizedImage;            
        end     
        
    end
    

end
