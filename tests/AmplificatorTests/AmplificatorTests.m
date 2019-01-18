classdef AmplificatorTests < testRunner
    
    
    properties (Access = protected)
        FieldOfStudy = 'Amplificator'
        tests
    end
    
    methods (Access = public)
        function obj = AmplificatorTests()
            obj@testRunner();
        end
    end
    
    methods (Access = protected)
        function loadTests(obj)
            obj.tests = {...   
                   'testAmplificatorTensorForInclusions';
                   'testAmplificatorTensorNumericVsExplicitForSeqLam';
                   'testHomogenizationLaminateForAnisotropicTensors';
                };

        end
    end
    
end

