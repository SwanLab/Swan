classdef testRunner < handle
    
    properties (Abstract, Access = protected)
        FieldOfStudy
        tests
    end
    
    methods (Access = protected)
        
        function obj = testRunner()
            obj.addPath()
            obj.printHeader()
            obj.compute()
            obj.printTail()
        end
        
        function runTests(obj)
            nTests = size(obj.tests,1);
            for itest = 1: nTests
                test = eval(obj.tests{itest});
                test.checkTestPassed(obj.tests{itest});
            end
        end
        
    end
    
    
    methods (Abstract, Access = protected)
        loadTests(obj)
    end
    
    methods (Access = private)
        
        function printHeader(obj)
            fprintf(['Running, ',obj.FieldOfStudy,' tests...\n'])
        end
        
        function printTail(obj)
            fprintf(['\n',obj.FieldOfStudy,' tests completed.\n'])
            fprintf('\n-------------------------------------------\n\n')
        end
        
        function compute(obj)
            obj.loadTests()
            obj.runTests()
        end
        
    end
    
    methods (Access = private, Static)
        
        function addPath()
            addpath(genpath(pwd));
        end
    end
    
end