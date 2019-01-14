classdef testRunner < handle
    properties (Abstract, Access = protected)
        FieldOfStudy
        tests
    end
    
    methods (Access = protected)
        
        function obj = testRunner()
            obj.printHeader()
            obj.compute()
            obj.printTail()
        end
        
        function compute(obj)
            obj.loadTests()
            obj.runTests()
        end
        
        function runTests(obj)
            for itest = 1: size(obj.tests,1)
                test = eval(obj.tests{itest});
                test.checkTestPassed(obj.tests{itest});
            end
        end
       
        function printHeader(obj)
            fprintf(['Running, ',obj.FieldOfStudy,' tests...\n'])
        end
        
        function printTail(obj)
            fprintf(['\n',obj.FieldOfStudy,' tests completed.\n'])
            fprintf('\n-------------------------------------------\n\n')
        end
    end
    
    
    methods (Abstract, Access = protected)
        loadTests(obj)
    end
    
    
end

