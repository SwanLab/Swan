classdef PerformanceTestsSuite < handle

    methods

        function obj = PerformanceTestsSuite()
            path = './tests/Source/PerformanceTests/PerformanceTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','Performance');
            results = suite.run;
            table(results)
        end

    end
end