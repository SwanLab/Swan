classdef testAmplificatorTensorForSquareInclusion < ...
        testShowingError & ...
        testLoadStoredVariable & ...
        testStoredComputedChecker
    
    
    properties (Access = private)
        ampTensorNum
        ampTensorExp
    end
    
    properties (Access = protected)
        testName = 'AmplificatorTensorForSquareInclusion';  
        variablesToStore = {'P'};
        tol = 1e-6; 
    end
       
    
    methods (Access = public)
        
        function obj = testAmplificatorTensorForSquareInclusion() 
            obj.createNumericalAmplificationTensor()
            obj.selectComputedVar()
        end
      
        
    end
    
    methods (Access = private)
        
                                
       function createNumericalAmplificationTensor(obj)
            printTopology  = false;
%             homogenizer    = NumericalSquareHomogenizer(printTopology);
%             Ch             = homogenizer.Ch;
%             P              = homogenizer.getAmplificatorTensor();
            obj.ampTensorNum = SymmetricFourthOrderPlaneStressVoigtTensor();
            obj.ampTensorNum.createRandomTensor();
        end
        

    end    
    
    methods (Access = protected)
        
        function selectComputedVar(obj)
            obj.computedVar{1} = obj.ampTensorNum.getValue();            
        end
        
    end
    

end

