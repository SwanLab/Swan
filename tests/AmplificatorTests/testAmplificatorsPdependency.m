classdef testAmplificatorsPdependency < ...
        testShowingError & ...
        testLoadStoredVariable & ...
        testStoredComputedChecker
    
    
    properties (Access = protected)
       variablesToStore = {'P'};
        tol = 1e-6;
        testName = 'AmplificatorsPdependency.mat';       
    end
    
    properties (Access = private)
       nExp
       pExp
       firstInvP
       m1 = 0.2;
       m2 = 0.8; 
       printing = false;
       iter = 0;
       hSmooth
       hNonSmooth
    end
    
    methods (Access = public)
        
        function obj = testAmplificatorsPdependency() 
            obj.init();
            obj.createSmoothRectInclusionHomogenizer();
            obj.createNonSmoothRectInclusionHomogenizer();
            obj.computeFirstAmplicatorComponents();
            obj.selectComputedVar();
            obj.plotFirstAmplificatorComponent();
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.nExp = 6;
            obj.pExp = 2.^[1:obj.nExp];            
        end
        
        function computeFirstAmplicatorComponents(obj)
            obj.firstInvP = ones(obj.nExp,2);
            obj.firstInvP(2,1) = obj.obtainFirstPinverseFromHomogenizer(obj.hSmooth);
            obj.firstInvP(2,2) = obj.obtainFirstPinverseFromHomogenizer(obj.hNonSmooth);
        end
        
        function selectComputedVar(obj)
            obj.computedVar{1} = obj.firstInvP;
        end
        
        function plotFirstAmplificatorComponent(obj)
            plot(obj.pExp,obj.firstInvP,'-+')
        end
        
        function createSmoothRectInclusionHomogenizer(obj)
            f    = strcat(obj.testName);
            p    = obj.printing;
            m1v  = obj.m1;
            m2v  = obj.m2; 
            i    = obj.iter;
            obj.hSmooth = NumericalSmoothRectangleHomogenizer(f,p,m1v,m2v,i);                
        end
        
        function createNonSmoothRectInclusionHomogenizer(obj)
            f    = strcat(obj.testName);
            p    = obj.printing;
            m1v  = obj.m1;
            m2v  = obj.m2; 
            i    = obj.iter;
            obj.hNonSmooth = NumericalRectangleHomogenizer(f,p,m1v,m2v,i);                
        end
        
        function f = obtainFirstPinverseFromHomogenizer(obj,homog)
            Ch            = homog.getCh();
            Ptensor       = homog.getAmplificatorTensor();
            PTensorValues = Ptensor.getValue();
            f = 1/PTensorValues(1,1);            
        end
        
    end
    
    
end