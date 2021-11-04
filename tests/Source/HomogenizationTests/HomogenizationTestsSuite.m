classdef HomogenizationTestsSuite < handle
    
    methods

        function obj = HomogenizationTestsSuite()
            path = './tests/Source/HomogenizationTests/HomogenizationTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','HomogenizationTests');
            results = suite.run;
            table(results)
        end

    end
end