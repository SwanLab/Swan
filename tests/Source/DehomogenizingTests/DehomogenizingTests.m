classdef DehomogenizingTests < testRunner
    
    properties (Access = protected)
        FieldOfStudy = 'Dehomogenizing'
        tests
    end
    
    methods (Access = public)
        
        function  obj = DehomogenizingTests()
            obj@testRunner();
        end
    end
    
    methods (Access = protected)
        function loadTests(obj)
            obj.tests = {'MeshSymmetrizerTest'; 
                'ScalarSymmetrizerTest'};            
        end
    end
end