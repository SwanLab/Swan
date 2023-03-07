classdef FunctionTestsSuite < handle
    
    methods

        function obj = FunctionTestsSuite()
            path = './tests/Source/FunctionTests/FunctionTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','FunctionTests');
            results = suite.run;
            table(results)
        end

    end
end