classdef testAmplificatorTensorForInclusions < ...
        testShowingError & ...
        testLoadStoredVariable & ...
        testStoredComputedChecker
    
    properties (Access = private)
        ampTensorNum
    end
    
    properties (Access = protected)
        variablesToStore = {'P'};
        tol = 1e-6;
        fileName
        printTopology        
    end
    
    methods (Access = protected)
        
        function obj = testAmplificatorTensorForInclusions()
            obj.createNumericalAmplificationTensor()
            obj.selectComputedVar()            
        end
        
        function selectComputedVar(obj)
            obj.computedVar{1} = obj.ampTensorNum.getValue();
        end
        
    end
    
    methods (Access = private)
        
        function createNumericalAmplificationTensor(obj)
            n = 5;
            m1vect = linspace(0.001,0.99,n);
            P = zeros(9,n);
            obj.printTopology  = true;
            obj.fileName       = strcat(obj.testName);
            for im1 = 1:n
                m1            = m1vect(im1);
                homog         = obj.createHomogenizer(im1,m1);
                Ch            = homog.getCh();
                Ptensor       = homog.getAmplificatorTensor();
                PTensorValues = Ptensor.getValue();
                P(:,im1)      = PTensorValues(:);
            end
            obj.ampTensorNum = SymmetricFourthOrderPlaneStressVoigtTensor();
            obj.ampTensorNum.createRandomTensor();
        end        
        
    end
    
    methods (Access = protected, Abstract)
        createHomogenizer(obj,iter,m1)        
    end
    
end

